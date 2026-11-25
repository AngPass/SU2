/*!
 * \file CFVMFlowSolverBase_LES.inl
 * \brief Finite volume LES closures (SOLVER= LES / INC_LES): explicit SGS model, stochastic backscatter
 *        (Leith) and LES wall model, with the wall model exchange location at the first node off the wall
 *        or at a height h_wm.
 * \note Included by CFVMFlowSolverBase.inl.
 * \version 8.6.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../../../Common/include/adt/CADTElemClass.hpp"
#include "../../../Common/include/toolboxes/random_toolbox.hpp"

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::InitLES(const CGeometry* geometry, const CConfig* config,
                                                unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- Explicit SGS model (validated in CConfig). The model constants are the ones of the FEM-DG solver,
   *    except for Smagorinsky whose filter multiplier is 1, i.e. C_s * Delta with C_s = 0.1 and Delta
   *    the finite volume filter width. ---*/

  if (config->GetFVM_SGS_Model()) {
    switch (config->GetKind_SGS_Model()) {
      case TURB_SGS_MODEL::SMAGORINSKY: {
        auto smag = std::make_unique<CSmagorinskyModel>();
        smag->filter_mult = 1.0;
        SGSModel = std::move(smag);
        break;
      }
      case TURB_SGS_MODEL::WALE:
        SGSModel = std::make_unique<CWALEModel>();
        break;
      case TURB_SGS_MODEL::VREMAN:
        SGSModel = std::make_unique<CVremanModel>();
        break;
      default:
        break;
    }
  }
  if (SGSModel) SGS_EddyViscosity.resize(nPoint) = su2double(0.0);

  /*--- Stochastic backscatter (Leith), validated in CConfig (WHITE_NOISE, explicit SGS model). ---*/

  if (config->GetLES_FVM() && config->GetSBSParam().StochasticBackscatter && iMesh == MESH_0) {
    if (nDim < 3)
      SU2_MPI::Error("The stochastic backscatter of SOLVER= LES / INC_LES is available in 3D only.", CURRENT_FUNCTION);
    nodes->AllocateLES_Backscatter();
  }

  /*--- LES wall model: modeled heat flux storage, exchange locations at h_wm > 0 (collective: the
   *    conditions are global). ---*/

  if (config->GetLES_FVM() && config->GetWall_Functions()) {
    WM_HeatFlux.resize(nMarker);
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) WM_HeatFlux[iMarker].assign(nVertex[iMarker], 0.0);

    if (iMesh == MESH_0) SetupWallModelExchange(geometry, config);
  }
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::ComputeSGS_EddyViscosity(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const su2double userFilterWidth = config->GetLES_FilterWidth();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Filter width: user value, or from the (complete, i.e. including periodic) dual volume. ---*/

    su2double filterWidth = userFilterWidth;
    if (filterWidth <= 0.0) {
      const su2double vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      filterWidth = (nDim == 3) ? pow(vol, 1.0/3.0) : sqrt(vol);
    }

    const su2double rho = nodes->GetDensity(iPoint);
    const su2double dist = geometry->nodes->GetWall_Distance(iPoint);
    const auto grad = nodes->GetVelocityGradient(iPoint);

    su2double muT = 0.0;
    if (nDim == 2) {
      muT = SGSModel->ComputeEddyViscosity_2D(rho, grad(0,0), grad(0,1), grad(1,0), grad(1,1),
                                              filterWidth, dist);
    } else {
      muT = SGSModel->ComputeEddyViscosity_3D(rho, grad(0,0), grad(0,1), grad(0,2),
                                              grad(1,0), grad(1,1), grad(1,2),
                                              grad(2,0), grad(2,1), grad(2,2),
                                              filterWidth, dist);
    }
    muT = max(muT, su2double(0.0));

    /*--- Stored for SetPrimitive_Variables (thermal conductivity, next iteration), and set
     *    directly in the primitives so that the current viscous residual uses it. ---*/

    SGS_EddyViscosity(iPoint) = muT;
    nodes->SetEddyViscosity(iPoint, muT);
  }
  END_SU2_OMP_FOR
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::SetLES_Backscatter(CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const unsigned long timeIter = config->GetTimeIter();
  const su2double dt = config->GetDelta_UnstTimeND();
  const unsigned short nPeriodic = config->GetnMarker_Periodic() / 2;
  const bool inBox = config->GetSBSParam().StochBackscatterInBox;
  const auto* boxBounds = config->GetSBSParam().StochBackscatterBoxBounds;

  auto PeriodicComm = [&](unsigned short commType) {
    for (unsigned short iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
      InitiatePeriodicComms(geometry, config, iPeriodic, commType);
      CompletePeriodicComms(geometry, config, iPeriodic, commType);
    }
  };

  /*--- 1) Unit Gaussian noise, white in space and time. The seed depends on the global index (independent of
   *       the partitioning); the passive periodic face takes the noise of the master face. ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    const unsigned long iPointGlobal = geometry->nodes->GetGlobalIndex(iPoint);
    for (unsigned short iDim = 0; iDim < 3; iDim++)
      nodes->SetSBSPeriodicBuf(iPoint, iDim, RandomToolbox::GetNormal(iPointGlobal, iDim, timeIter));
  }
  END_SU2_OMP_FOR

  PeriodicComm(PERIODIC_SBS_COPY);

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (unsigned short iDim = 0; iDim < 3; iDim++)
      nodes->SetLangevinSourceTerms(iPoint, iDim, nodes->GetSBSPeriodicBuf(iPoint, iDim));
  END_SU2_OMP_FOR

  InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
  CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);

  /*--- 2) Smoothing filter over the nearest neighbours (Leith's 2 dx filter): xi_i + sum_j xi_j, partial sums
   *       completed at periodic points. ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double sum[3] = {0.0, 0.0, 0.0};
    const auto nNeigh = geometry->nodes->GetnPoint(iPoint);
    for (unsigned short iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
      const auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
      for (unsigned short iDim = 0; iDim < 3; iDim++) sum[iDim] += nodes->GetLangevinSourceTerms(jPoint, iDim);
    }
    for (unsigned short iDim = 0; iDim < 3; iDim++) nodes->SetSBSPeriodicBuf(iPoint, iDim, sum[iDim]);
    nodes->SetSBSPeriodicBuf(iPoint, 3, su2double(nNeigh));
  }
  END_SU2_OMP_FOR

  PeriodicComm(PERIODIC_SBS_SUM);

  /*--- 3) Vector potential phi = nu_sgs sqrt(|S|/dt) xi, with xi = (xi_i + sum_j xi_j) / sqrt(1 + n), of unit
   *       variance (exactly, away from periodic faces, the samples being independent). This is Leith's
   *       C_b |S dt|^(3/2) (lambda/dt)^2 with the Smagorinsky estimate (C_s lambda)^2 |S| replaced by the SGS eddy
   *       viscosity, so that the forcing vanishes with nu_sgs at the walls (WALE, Vreman). ---*/

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    const auto grad = nodes->GetVelocityGradient(iPoint);
    su2double div = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) div += grad(iDim, iDim);
    su2double S2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short jDim = 0; jDim < nDim; jDim++) {
        const su2double s_ij = 0.5 * (grad(iDim, jDim) + grad(jDim, iDim)) - ((iDim == jDim) ? div / 3.0 : 0.0);
        S2 += 2.0 * s_ij * s_ij;
      }
    }
    const su2double nuSGS = SGS_EddyViscosity(iPoint) / nodes->GetDensity(iPoint);
    su2double amplitude = nuSGS * sqrt(sqrt(S2) / dt);

    /*--- SBS_IN_BOX: sharp mask on the smoothed field (unit variance up to the box faces), evaluated at every
     *    time step to follow a moving grid. ---*/

    if (inBox) {
      const auto coord = geometry->nodes->GetCoord(iPoint);
      for (unsigned short iDim = 0; iDim < 3; iDim++)
        if (coord[iDim] < boxBounds[2*iDim] || coord[iDim] > boxBounds[2*iDim+1]) amplitude = 0.0;
    }
    const su2double norm = 1.0 / sqrt(1.0 + nodes->GetSBSPeriodicBuf(iPoint, 3));

    for (unsigned short iDim = 0; iDim < 3; iDim++) {
      const su2double xi = (nodes->GetLangevinSourceTerms(iPoint, iDim) + nodes->GetSBSPeriodicBuf(iPoint, iDim)) * norm;
      nodes->SetLangevinSourceTerms(iPoint, iDim, amplitude * xi);
    }
  }
  END_SU2_OMP_FOR

  /*--- 4) No forcing on the physical boundaries. ---*/

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    const auto kindBC = config->GetMarker_All_KindBC(iMarker);
    if (kindBC == SEND_RECEIVE || kindBC == PERIODIC_BOUNDARY || kindBC == INTERNAL_BOUNDARY ||
        kindBC == NEARFIELD_BOUNDARY || kindBC == FLUID_INTERFACE) continue;

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      for (unsigned short iDim = 0; iDim < 3; iDim++) nodes->SetLangevinSourceTerms(iPoint, iDim, 0.0);
    }
    END_SU2_OMP_FOR
  }

  /*--- Halo values, read by the viscous residual. ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
  CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::SetupWallModelExchange(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const auto nMarkerLocal = config->GetnMarker_All();
  const int size = SU2_MPI::GetSize();
  const int rank = SU2_MPI::GetRank();
  const auto nVarExch = WM_nVar();

  /*--- Exchange height of each LES wall model marker; is it used anywhere (on any rank)? ---*/

  WM_Height.assign(nMarkerLocal, 0.0);
  int localHeight = 0, globalHeight = 0;
  for (auto iMarker = 0u; iMarker < nMarkerLocal; iMarker++) {
    if (!config->GetViscous_Wall(iMarker)) continue;
    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::LOGARITHMIC_MODEL) continue;
    WM_Height[iMarker] = config->GetWallFunction_DoubleInfo(Marker_Tag)[0];
    if (WM_Height[iMarker] > 0.0) localHeight = 1;
  }
  SU2_MPI::Allreduce(&localHeight, &globalHeight, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
  if (globalHeight == 0) return;

  if (dynamic_grid)
    SU2_MPI::Error("LES wall model with h_wm > 0 is not available with grid motion.", CURRENT_FUNCTION);

  WM_Exchange = true;

  /*--- Local ADT of the primal elements, including the halo ones (their nodes are all local). ---*/

  const auto nPointLocal = geometry->GetnPoint();
  const auto nElemLocal = geometry->GetnElem();

  vector<su2double> coor(nPointLocal * nDim);
  for (auto iPoint = 0ul; iPoint < nPointLocal; iPoint++)
    for (auto iDim = 0u; iDim < nDim; iDim++) coor[iPoint*nDim + iDim] = geometry->nodes->GetCoord(iPoint, iDim);

  vector<unsigned long> elemConn, elemIDs(nElemLocal);
  vector<unsigned short> elemVTK(nElemLocal), elemMarker(nElemLocal, 0);
  for (auto iElem = 0ul; iElem < nElemLocal; iElem++) {
    const auto elem = geometry->elem[iElem];
    elemVTK[iElem] = elem->GetVTK_Type();
    elemIDs[iElem] = iElem;
    for (auto iNode = 0u; iNode < elem->GetnNodes(); iNode++) elemConn.push_back(elem->GetNode(iNode));
  }

  CADTElemClass localADT(nDim, coor, elemConn, elemVTK, elemMarker, elemIDs, false);

  auto search = [&](const su2double* x, CWMDonor& donor) {
    unsigned short markerID;
    unsigned long elemID;
    int rankID;
    su2double parCoor[3] = {0.0}, weights[8] = {0.0};
    if (!localADT.DetermineContainingElement(x, markerID, elemID, rankID, parCoor, weights)) return false;
    const auto elem = geometry->elem[elemID];
    donor.nNodes = elem->GetnNodes();
    for (auto iNode = 0u; iNode < donor.nNodes; iNode++) {
      donor.node[iNode] = elem->GetNode(iNode);
      donor.weight[iNode] = weights[iNode];
    }
    return true;
  };

  /*--- Exchange locations x_e = x_w + h_wm * n (n: inward unit normal) of all the wall vertices, halo ones
   *    included (as for the first-node model, Tau_Wall is then available on halos without communication). ---*/

  WM_Index.resize(nMarkerLocal);
  vector<su2double> remoteCoor;   // exchange locations without a local donor
  vector<unsigned long> remoteK;  // their exchange location index
  unsigned long nExchange = 0;

  for (auto iMarker = 0u; iMarker < nMarkerLocal; iMarker++) {
    WM_Index[iMarker].assign(geometry->nVertex[iMarker], -1);
    if (!(WM_Height[iMarker] > 0.0)) continue;

    for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      const su2double Area = GeometryToolbox::Norm(nDim, Normal);

      su2double xe[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        xe[iDim] = geometry->nodes->GetCoord(iPoint, iDim) - WM_Height[iMarker] * Normal[iDim] / Area;

      WM_Index[iMarker][iVertex] = nExchange;
      CWMDonor donor;
      if (!search(xe, donor)) {
        remoteCoor.insert(remoteCoor.end(), xe, xe + nDim);
        remoteK.push_back(nExchange);
      }
      WM_LocalDonor.push_back(donor);
      nExchange++;
    }
  }

  WM_ExchangeVar.assign(nExchange * nVarExch, 0.0);
  WM_RecvSlot.assign(nExchange, -1);
  WM_SendCounts.assign(size, 0); WM_SendDispls.assign(size, 0);
  WM_RecvCounts.assign(size, 0); WM_RecvDispls.assign(size, 0);

  /*--- Exchange locations outside the local partition: every rank searches the locations of the others,
   *    the lowest rank containing a location serves it. ---*/

  const int nRemote = remoteK.size();
  vector<int> nRemoteAll(size, 0);
  SU2_MPI::Allgather(&nRemote, 1, MPI_INT, nRemoteAll.data(), 1, MPI_INT, SU2_MPI::GetComm());

  vector<int> ptOffset(size + 1, 0), coorCounts(size), coorDispls(size);
  for (int r = 0; r < size; r++) {
    ptOffset[r + 1] = ptOffset[r] + nRemoteAll[r];
    coorCounts[r] = nRemoteAll[r] * nDim;
    coorDispls[r] = ptOffset[r] * nDim;
  }
  const int nRemoteTotal = ptOffset[size];

  if (nRemoteTotal > 0) {
    vector<su2double> allCoor(nRemoteTotal * nDim);
    SU2_MPI::Allgatherv(remoteCoor.data(), nRemote * nDim, MPI_DOUBLE, allCoor.data(), coorCounts.data(),
                        coorDispls.data(), MPI_DOUBLE, SU2_MPI::GetComm());

    /*--- Search the locations of the other ranks (in rank and then location order). ---*/

    struct CCandidate { int global; int requester; CWMDonor donor; };
    vector<CCandidate> candidates;
    vector<int> foundRank(nRemoteTotal, size), owner(nRemoteTotal, size);

    for (int r = 0; r < size; r++) {
      if (r == rank) continue;
      for (int g = ptOffset[r]; g < ptOffset[r + 1]; g++) {
        CWMDonor donor;
        if (search(&allCoor[g * nDim], donor)) {
          foundRank[g] = rank;
          candidates.push_back({g, r, donor});
        }
      }
    }
    SU2_MPI::Allreduce(foundRank.data(), owner.data(), nRemoteTotal, MPI_INT, MPI_MIN, SU2_MPI::GetComm());

    /*--- Serving side: the candidates this rank owns, grouped by requesting rank in location order. ---*/

    for (const auto& c : candidates) {
      if (owner[c.global] != rank) continue;
      WM_Serve.push_back(c.donor);
      WM_SendCounts[c.requester]++;
    }

    /*--- Requesting side: receive slot of each remote location, grouped by owner in location order. ---*/

    for (int j = 0; j < nRemote; j++) {
      const int o = owner[ptOffset[rank] + j];
      if (o == size)
        SU2_MPI::Error("LES wall model: an exchange location at h_wm from the wall is outside the domain, "
                       "reduce h_wm.", CURRENT_FUNCTION);
      WM_RecvCounts[o]++;
    }
    for (int r = 1; r < size; r++) {
      WM_SendDispls[r] = WM_SendDispls[r - 1] + WM_SendCounts[r - 1];
      WM_RecvDispls[r] = WM_RecvDispls[r - 1] + WM_RecvCounts[r - 1];
    }
    vector<int> cursor(size, 0);
    for (int j = 0; j < nRemote; j++) {
      const int o = owner[ptOffset[rank] + j];
      WM_RecvSlot[remoteK[j]] = WM_RecvDispls[o] + cursor[o]++;
    }

    /*--- From locations to values for Alltoallv. ---*/

    for (int r = 0; r < size; r++) {
      WM_SendCounts[r] *= nVarExch; WM_SendDispls[r] *= nVarExch;
      WM_RecvCounts[r] *= nVarExch; WM_RecvDispls[r] *= nVarExch;
    }
  }

  /*--- Summary. ---*/

  unsigned long nLocal[2] = {nExchange, static_cast<unsigned long>(nRemote)}, nGlobal[2] = {0, 0};
  SU2_MPI::Allreduce(nLocal, nGlobal, 2, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  if (rank == MASTER_NODE)
    cout << "LES wall model: " << nGlobal[0] << " exchange locations at h_wm (halo vertices included), "
         << nGlobal[1] << " served by another rank." << endl;
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::UpdateWallModelExchange() {
  SU2_ZONE_SCOPED

  const auto nVarExch = WM_nVar();

  /*--- Velocity components, then temperature. ---*/

  auto interpolate = [&](const CWMDonor& donor, su2double* var) {
    for (auto iVar = 0u; iVar < nVarExch; iVar++) var[iVar] = 0.0;
    for (auto iNode = 0u; iNode < donor.nNodes; iNode++) {
      const auto jPoint = donor.node[iNode];
      for (auto iDim = 0u; iDim < nDim; iDim++) var[iDim] += donor.weight[iNode] * nodes->GetVelocity(jPoint, iDim);
      var[nDim] += donor.weight[iNode] * nodes->GetTemperature(jPoint);
    }
  };

  /*--- Local donors. ---*/

  for (auto k = 0ul; k < WM_LocalDonor.size(); k++)
    if (WM_LocalDonor[k].nNodes > 0) interpolate(WM_LocalDonor[k], &WM_ExchangeVar[k * nVarExch]);

  /*--- Remote donors (collective, all ranks take part even with nothing to send or receive). ---*/

  const int size = SU2_MPI::GetSize();
  if (size == 1) return;

  const int nSend = WM_SendDispls[size - 1] + WM_SendCounts[size - 1];
  const int nRecv = WM_RecvDispls[size - 1] + WM_RecvCounts[size - 1];

  /*--- Buffers of at least one entry, the call is collective even for ranks with nothing to exchange. ---*/
  vector<su2double> sendBuf(max(nSend, 1)), recvBuf(max(nRecv, 1));
  for (auto i = 0ul; i < WM_Serve.size(); i++) interpolate(WM_Serve[i], &sendBuf[i * nVarExch]);

  SU2_MPI::Alltoallv(sendBuf.data(), WM_SendCounts.data(), WM_SendDispls.data(), MPI_DOUBLE, recvBuf.data(),
                     WM_RecvCounts.data(), WM_RecvDispls.data(), MPI_DOUBLE, SU2_MPI::GetComm());

  for (auto k = 0ul; k < WM_RecvSlot.size(); k++) {
    if (WM_RecvSlot[k] < 0) continue;
    for (auto iVar = 0u; iVar < nVarExch; iVar++)
      WM_ExchangeVar[k * nVarExch + iVar] = recvBuf[WM_RecvSlot[k] * nVarExch + iVar];
  }
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::SetTau_Wall_LES(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- Reichardt's law of the wall, valid from the viscous sublayer to the log layer:
   *      u+ = ln(1 + kappa y+)/kappa + A (1 - exp(-y+/11) - (y+/11) exp(-0.33 y+)),  A = B - ln(kappa)/kappa,
   *    with B the log-law constant (WALLMODEL_B, u+ -> ln(y+)/kappa + B). Since it also covers y+ -> 0
   *    (u+ -> y+), no switch-off is needed when the exchange location falls in the viscous sublayer. ---*/

  const su2double kappa = config->GetwallModel_Kappa();
  const su2double A = config->GetwallModel_B() - log(kappa) / kappa;
  const auto max_iter = config->GetwallModel_MaxIter();
  const su2double Prandtl_Lam = config->GetPrandtl_Lam();
  const auto nVarExch = WM_nVar();
  constexpr passivedouble tol = 1e-10;

  auto uPlus = [&](su2double yp) {
    return log(1.0 + kappa * yp) / kappa + A * (1.0 - exp(-yp / 11.0) - (yp / 11.0) * exp(-0.33 * yp));
  };
  auto duPlus_dyPlus = [&](su2double yp) {
    return 1.0 / (1.0 + kappa * yp) +
           A * (exp(-yp / 11.0) / 11.0 - exp(-0.33 * yp) / 11.0 + 0.33 * (yp / 11.0) * exp(-0.33 * yp));
  };

  /*--- Kader's law for the temperature, T+ = (T_w - T) rho_w cp u_tau / q_w:
   *      T+ = Pr y+ exp(G) + (2.12 ln(1 + y+) + beta) exp(1/G),
   *      G = -0.01 (Pr y+)^4 / (1 + 5 Pr^3 y+),  beta = (3.85 Pr^(1/3) - 1.3)^2 + 2.12 ln(Pr). ---*/

  auto tPlus = [&](su2double yp) {
    const su2double G = -0.01 * pow(Prandtl_Lam * yp, 4.0) / (1.0 + 5.0 * pow(Prandtl_Lam, 3.0) * yp);
    const su2double beta = pow(3.85 * pow(Prandtl_Lam, 1.0/3.0) - 1.3, 2.0) + 2.12 * log(Prandtl_Lam);
    return Prandtl_Lam * yp * exp(G) + (2.12 * log(1.0 + yp) + beta) * exp(1.0 / fmin(G, -EPS));
  };

  /*--- Solution at the exchange locations at h_wm (interpolation + MPI, one thread, all ranks). ---*/

  const bool exchangeAtHeight = WM_Exchange && (MGLevel == MESH_0);
  if (exchangeAtHeight) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    UpdateWallModelExchange();
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    if (!config->GetViscous_Wall(iMarker)) continue;

    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::LOGARITHMIC_MODEL) continue;

    /*--- Modeled heat flux on isothermal walls of the compressible solver (the incompressible one imposes
     *    the wall temperature strongly, its wall heat flux follows from the resolved gradient). ---*/

    const bool heatFlux = (FlowRegime == ENUM_REGIME::COMPRESSIBLE) &&
                          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) &&
                          !config->GetMarker_All_PyCustom(iMarker);
    const su2double Twall = heatFlux ? config->GetIsothermal_Temperature(Marker_Tag) / config->GetTemperature_Ref()
                                     : su2double(0.0);

    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

      /*--- On the finest mesh compute also on halo nodes to avoid communication of tau wall. ---*/
      if ((!geometry->nodes->GetDomain(iPoint)) && !(MGLevel==MESH_0)) continue;

      /*--- Inward unit normal of the wall. ---*/

      const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      const su2double Area = GeometryToolbox::Norm(nDim, Normal);
      su2double UnitNormal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Exchange location and instantaneous solution there (velocity relative to the wall):
       *    - h_wm > 0: location at distance h_wm along the normal, interpolated in its donor element;
       *    - h_wm <= 0: first node off the wall, distance as the normal projection max(|n.d|, 0.05 |d|). ---*/

      const long kExchange = exchangeAtHeight ? WM_Index[iMarker][iVertex] : -1;

      su2double y = 0.0, Texch = 0.0;
      su2double Vel[MAXNDIM] = {0.0};

      if (kExchange >= 0) {
        y = WM_Height[iMarker];
        for (auto iDim = 0u; iDim < nDim; iDim++) Vel[iDim] = WM_ExchangeVar[kExchange * nVarExch + iDim];
        Texch = WM_ExchangeVar[kExchange * nVarExch + nDim];
      } else {
        y = GeometryToolbox::NormalDistance(nDim, UnitNormal, geometry->nodes->GetCoord(Point_Normal),
                                            geometry->nodes->GetCoord(iPoint));
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          Vel[iDim] = nodes->GetVelocity(Point_Normal, iDim);
          if (dynamic_grid) Vel[iDim] -= geometry->nodes->GetGridVel()(iPoint, iDim);
        }
        Texch = nodes->GetTemperature(Point_Normal);
      }
      const su2double VelNormal = GeometryToolbox::DotProduct(nDim, Vel, UnitNormal);
      for (auto iDim = 0u; iDim < nDim; iDim++) Vel[iDim] -= VelNormal * UnitNormal[iDim];
      const su2double VelTang = GeometryToolbox::Norm(nDim, Vel);

      /*--- Wall fluid properties (wall node). ---*/

      const su2double rho = nodes->GetDensity(iPoint);
      const su2double nu = nodes->GetLaminarViscosity(iPoint) / rho;

      /*--- Newton iterations on g(u_tau) = u_tau u+(y u_tau / nu) - U = 0, monotonically increasing.
       *    Initial guess: previous value, or the laminar (viscous sublayer) estimate. ---*/

      su2double uTau = UTau[iMarker][iVertex];
      if (!(uTau > 0.0)) uTau = sqrt(nu * VelTang / y);
      uTau = max(uTau, su2double(1e-12));

      for (auto iter = 0u; iter < max_iter; iter++) {
        const su2double yPlus = y * uTau / nu;
        const su2double g = uTau * uPlus(yPlus) - VelTang;
        const su2double dg = uPlus(yPlus) + yPlus * duPlus_dyPlus(yPlus);
        su2double uTauNew = uTau - g / dg;
        if (uTauNew <= 0.0) uTauNew = 0.5 * uTau;  // keep positive
        const su2double change = fabs(uTauNew - uTau);
        uTau = uTauNew;
        if (change <= tol * uTau) break;
      }

      /*--- Store. Tau_Wall > 0 marks the node for the viscous flux, which then imposes this value on the
       *    wall-normal edges (zero tangential velocity gives a vanishing, but still imposed, stress). ---*/

      const su2double yPlus = y * uTau / nu;
      UTau[iMarker][iVertex] = uTau;
      YPlus[iMarker][iVertex] = yPlus;
      nodes->SetTau_Wall(iPoint, rho * uTau * uTau);

      /*--- Wall heat flux into the fluid, q_w = (T_w - T) rho cp u_tau / T+. ---*/

      if (heatFlux) {
        const su2double cp = nodes->GetSpecificHeatCp(iPoint);
        WM_HeatFlux[iMarker][iVertex] = (Twall - Texch) * rho * cp * uTau / tPlus(yPlus);
      }
    }
    END_SU2_OMP_FOR
  }
}
