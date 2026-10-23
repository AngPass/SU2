/*!
 * \file CTurbSASolver.cpp
 * \brief Main subroutines of CTurbSASolver class
 * \author F. Palacios, A. Bueno
 * \version 8.5.0 "Harrier"
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

#include "../../include/solvers/CTurbSASolver.hpp"
#include "../../include/solvers/CScalarSolver.inl"
#include "../../include/variables/CTurbSAVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../include/numerics/turbulent/turb_sa_edge_flux.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/random_toolbox.hpp"


CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, const CSolver* flow_solver, unsigned short iMesh,
                             CFluidModel* FluidModel)
             : CTurbSolver(geometry, config, flow_solver, false) {
  SU2_ZONE_SCOPED
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf, Factor_nu_Inf, Factor_nu_Engine, Factor_nu_ActDisk;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> dependent of the turbulent model ---*/

  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Add Langevin equations if the Stochastic Backscatter Model is used ---*/

  if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
    nVar += 3;
    nVarGrad = nPrimVar = nVar;

    /*--- Only nu_tilde (index 0) is MUSCL-reconstructed when MUSCL_TURB=YES; the 3 Langevin
          components always keep their raw nodal values, see CScalarFlux_SA::NVarReconstructed. ---*/
  }

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == MG_CYCLE::FULL) {

    /*--- Define some auxiliar vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy) {
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
      EdgeFluxesDiff.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
    }

    if (config->GetExtraOutput()) {
      if (nDim == 2) { nOutputVariables = 13; }
      else if (nDim == 3) { nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Read farfield conditions from config ---*/

  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();

  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/

  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  su2double nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  if (config->GetSAParsedOptions().bc) {
    nu_tilde_Inf  = 0.005*Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  }

  Solution_Inf[0] = nu_tilde_Inf;
  if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
    for (unsigned short iVar = 1; iVar < nVar; iVar++) {
      Solution_Inf[iVar] = 0.0;
    }
  }

  /*--- Factor_nu_Engine ---*/
  Factor_nu_Engine   = config->GetNuFactor_Engine();
  nu_tilde_Engine[0] = Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  if (config->GetSAParsedOptions().bc) {
    nu_tilde_Engine[0] = 0.005*Factor_nu_Engine*Viscosity_Inf/Density_Inf;
  }

  /*--- Factor_nu_ActDisk ---*/
  Factor_nu_ActDisk   = config->GetNuFactor_Engine();
  nu_tilde_ActDisk[0] = Factor_nu_ActDisk*Viscosity_Inf/Density_Inf;

  /*--- Eddy viscosity at infinity ---*/
  su2double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
  su2double muT_Inf;
  Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
  Ji_3 = Ji*Ji*Ji;
  fv1 = Ji_3/(Ji_3+cv1_3);
  muT_Inf = Density_Inf*fv1*nu_tilde_Inf;

  if (config->GetSAParsedOptions().version != SA_OPTIONS::NEG) {
    lowerlimit[0] = EPS;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Ghost states for boundary conditions, sized to the largest marker (see BoundaryFluxResidual). ---*/
  unsigned long maxMarkerVertices = 0;
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    maxMarkerVertices = max(maxMarkerVertices, nVertex[iMarker]);
  ghostNodes = make_unique<CTurbSAVariable>(nu_tilde_Inf, muT_Inf, maxMarkerVertices, nDim, nVar, config, true);

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar) = nu_tilde_Inf;
    if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
      for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        for (unsigned short iVar = 1; iVar < nVar; iVar++) {
          Inlet_TurbVars[iMarker](iVertex,iVar) = 0.0;
        }
      }
    }
  }

  /*--- Store the initial CFL number for all grid points. ---*/

  const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  CheckSBSSetup(config);

  /*--- Add the solver name. ---*/
  SolverName = "SA";

}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
        unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  const auto kind_hybridRANSLES = config->GetKind_HybridRANSLES();

  /*--- Clear Residual and Jacobian. Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);

  /*--- Undivided Laplacian of the 4th order artificial dissipation of the Langevin equations. ---*/
  const auto& sbs = config->GetSBSParam();
  if (!Output && sbs.StochasticBackscatter && sbs.stochSourceType == LANGEVIN && sbs.SBS_Kappa4 > 0.0)
    SetUndivided_Laplacian(geometry, config);

  if (kind_hybridRANSLES != NO_HYBRIDRANSLES) {

    /*--- Set the vortex tilting coefficient at every node if required ---*/

    if (kind_hybridRANSLES == SA_EDDES){
      auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        auto Vorticity = flowNodes->GetVorticity(iPoint);
        auto PrimGrad_Flow = flowNodes->GetGradient_Primitive(iPoint);
        auto Laminar_Viscosity = flowNodes->GetLaminarViscosity(iPoint);
        nodes->SetVortex_Tilting(iPoint, PrimGrad_Flow, Vorticity, Laminar_Viscosity);
      }
      END_SU2_OMP_FOR
    }

    SetDES_LengthScale(solver_container, geometry, config);

    bool backscatter = config->GetSBSParam().StochasticBackscatter;
    bool backscatterInBox = config->GetSBSParam().StochBackscatterInBox;

    /*--- maxDelta must reach halo points for SetStochSourceMom. ---*/
    if (backscatter) {
      InitiateComms(geometry, config, MPI_QUANTITIES::DES_FILTERWIDTH);
      CompleteComms(geometry, config, MPI_QUANTITIES::DES_FILTERWIDTH);
    }

    /*--- Only needed by backscatter/ hybrid models with filtered stresses. ---*/
    if (backscatter || config->GetSBSParam().filterStresses) {
      InitiateComms(geometry, config, MPI_QUANTITIES::LES_SENSOR);
      CompleteComms(geometry, config, MPI_QUANTITIES::LES_SENSOR);
    }

    /*--- Compute source terms for Langevin equations ---*/

    unsigned long innerIter = config->GetInnerIter();
    if (backscatter && innerIter==0) {
      /*--- Only read by SetLangevinSourceTerms, once per physical time step (kept here rather than
            in the constructor so that it follows the grid if it moves). ---*/
      if (backscatterInBox) SetBackscatterInBox(config, geometry);

      if (config->GetSBSParam().useMeanTurb) {
        InitiateComms(geometry, config, MPI_QUANTITIES::MEAN_EDDY_VISC);
        CompleteComms(geometry, config, MPI_QUANTITIES::MEAN_EDDY_VISC);
      }

      SetLangevinSourceTerms(config, geometry);
      const unsigned short maxIter = config->GetSBSParam().SBS_maxIterSmooth;
      if (maxIter > 0) SmoothLangevinSourceTerms(config, geometry);

      /*--- The final (smoothed and rescaled) source is computed on domain points only, but with
            WHITE_NOISE it is read at halo points too, by the flow viscous residual. ---*/
      InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
      CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);

      if (config->GetSBSParam().stochSourceType == ORNSTEIN_UHLENBECK) ComputeOU_Process(solver_container, config, geometry);
    }

  }

}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const su2double cv1_3 = 7.1*7.1*7.1, cR1 = 0.5, rough_const = 0.03;

  const bool neg_spalart_allmaras = config->GetSAParsedOptions().version == SA_OPTIONS::NEG;
  const bool backscatter = config->GetSBSParam().StochasticBackscatter;

  /*--- Compute eddy viscosity ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    const su2double rho = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    const su2double mu = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    const su2double nu = mu/rho;
    const su2double nu_hat = nodes->GetSolution(iPoint,0);
    const su2double roughness = geometry->nodes->GetRoughnessHeight(iPoint);
    const su2double dist = geometry->nodes->GetWall_Distance(iPoint) + rough_const * roughness;

    su2double Ji = nu_hat/nu;
    if (roughness > 1.0e-10)
      Ji += cR1*roughness/(dist+EPS);

    const su2double Ji_3 = Ji*Ji*Ji;
    const su2double fv1 = Ji_3/(Ji_3+cv1_3);

    su2double muT = rho*fv1*nu_hat;

    if ((neg_spalart_allmaras || backscatter) && nu_hat < 0) muT = 0.0;

    nodes->SetmuT(iPoint,muT);

  }
  END_SU2_OMP_FOR


  /*--- Compute turbulence index ---*/
  if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE || config->GetSAParsedOptions().bc) {
    auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

    for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      if (config->GetViscous_Wall(iMarker)) {
        SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
        for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

          if (geometry->nodes->GetDomain(iPoint)) {
            const auto jPoint = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            su2double FrictionVelocity = 0.0;
            /*--- Formulation varies for 2D and 3D problems: in 3D the friction velocity is assumed to be sqrt(mu * |Omega|)
            (provided by the reference paper https://doi.org/10.2514/6.1992-439), whereas in 2D we have to use the
            standard definition sqrt(c_f / rho) since Omega = 0.  ---*/
            if(nDim == 2){
              su2double shearStress = 0.0;
              for(auto iDim = 0u; iDim < nDim; iDim++) {
                shearStress += pow(solver_container[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, iDim), 2.0);
              }
              shearStress = sqrt(shearStress);

              FrictionVelocity = sqrt(shearStress/flowNodes->GetDensity(iPoint));
            } else {
              su2double VorticityMag = max(GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint)), 1e-12);
              FrictionVelocity = sqrt(flowNodes->GetLaminarViscosity(iPoint)*VorticityMag);
            }

            const su2double wall_dist = geometry->nodes->GetWall_Distance(jPoint);
            const su2double Derivative = nodes->GetSolution(jPoint, 0) / wall_dist;
            const su2double turbulence_index = Derivative / (FrictionVelocity * 0.41);

            nodes->SetTurbIndex(iPoint, turbulence_index);

          }
        }
        END_SU2_OMP_FOR
      }
    }
  }

  AD::EndNoSharedReading();
}

void CTurbSASolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics**,
                                    CConfig* config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const auto opt = ScalarFluxOptions::Interior(*config, config->GetBounded_Turb(),
                                               config->GetUse_Accurate_Turb_Jacobians());

  /*--- nVar is 1, or 4 with the three Langevin equations of stochastic backscatter. ---*/
  DispatchScheme<CScalarFlux_SA, 1, 4>(config, [&](auto tag) {
    EdgeFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt);
  });
}

void CTurbSASolver::BoundaryFlux(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                 const ScalarFluxOptions& opt, unsigned short val_marker) {
  DispatchScheme<CScalarFlux_SA, 1, 4>(config, [&](auto tag) {
    BoundaryFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, opt, val_marker);
  });
}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics*, CNumerics*,
                                 CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, Solution_Inf[iVar]);

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));
    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  bool axisymmetric = config->GetAxisymmetric();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  const bool transition_BC = config->GetSAParsedOptions().bc;

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  AD::StartNoSharedReading();

  /*--- Loop over all points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Get Hybrid RANS/LES Type and set the appropriate wall distance ---*/

    if (config->GetKind_HybridRANSLES() == NO_HYBRIDRANSLES) {

    /*--- For the SA model, wall roughness is accounted by modifying the computed wall distance
       *                              d_new = d + 0.03 k_s
       *    where k_s is the equivalent sand grain roughness height that is specified in cfg file.
       *    For smooth walls, wall roughness is zero and computed wall distance remains the same. */

      su2double modifiedWallDistance = geometry->nodes->GetWall_Distance(iPoint);

      modifiedWallDistance += 0.03*geometry->nodes->GetRoughnessHeight(iPoint);

      /*--- Set distance to the surface ---*/

      numerics->SetDistance(modifiedWallDistance, 0.0);

      /*--- Set the roughness of the closest wall. ---*/

      numerics->SetRoughness(geometry->nodes->GetRoughnessHeight(iPoint), 0.0 );

    } else {

      /*--- Set DES length scale ---*/

      numerics->SetDistance(nodes->GetDES_LengthScale(iPoint), 0.0);

      /*--- Compute source terms in Langevin equations (Stochastic Basckscatter Model) ---*/

      if (config->GetSBSParam().StochasticBackscatter) {
        if (config->GetSBSParam().stochSourceType == LANGEVIN || config->GetSBSParam().stochSourceType == WHITE_NOISE) {
          for (unsigned short iDim = 0; iDim < nDim; iDim++)
            numerics->SetStochSource(nodes->GetLangevinSourceTerms(iPoint, iDim), iDim);
        } else {
          for (unsigned short iDim = 0; iDim < nDim; iDim++)
            numerics->SetStochSource(nodes->GetOU_Process(iPoint, iDim), iDim);
        }
        numerics->SetLES_Mode(nodes->GetLES_Mode(iPoint), 0.0);
        numerics->SetMaxDelta(nodes->GetDES_FilterWidth(iPoint), 0.0);
        numerics->SetWallDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

        /*--- Time-averaged eddy viscosity, used instead of the instantaneous one to scale the
              stochastic forcing when SBS_USE_MEAN_TURB is active. ---*/
        if (config->GetSBSParam().useMeanTurb)
          numerics->SetAvgEddyViscosity(nodes->GetMeanEddyViscosity(iPoint), 0.0);
      }

    }

    /*--- Effective Intermittency ---*/

    if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
      numerics->SetIntermittencyEff(solver_container[TRANS_SOL]->GetNodes()->GetIntermittencyEff(iPoint));
      numerics->SetIntermittency(solver_container[TRANS_SOL]->GetNodes()->GetSolution(iPoint, 0));
    }

    if (axisymmetric) {
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));
    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Store the intermittency ---*/

    if (transition_BC || config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
      nodes->SetIntermittency(iPoint,numerics->GetIntermittencyEff());
    }
    
    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);

    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  if (harmonic_balance) {

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double Volume = geometry->nodes->GetVolume(iPoint);

      /*--- Access stored harmonic balance source term ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double Source = nodes->GetHarmonicBalance_Source(iPoint,iVar);
        LinSysRes(iPoint,iVar) += Source*Volume;
      }
    }
    END_SU2_OMP_FOR
  }

  AD::EndNoSharedReading();

  /*--- Custom user defined source term (from the python wrapper) ---*/
  if (config->GetPyCustomSource()) {
    CustomSourceResidual(geometry, solver_container, numerics_container, config, iMesh);
  }

}

void CTurbSASolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  /*--- Evaluate nu tilde at the closest point to the surface using the wall functions. ---*/

  if (config->GetWall_Functions()) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(SetTurbVars_WF(geometry, solver_container, config, val_marker);)
    return;
  }

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  WALL_TYPE WallType;
  su2double Roughness_Height;
  tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
  Roughness_Height = max(Roughness_Height, EPS);

  /*--- The dirichlet condition is used only without wall function, otherwise the
   convergence is compromised as we are providing nu tilde values for the
   first point of the wall  ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    if (WallType == WALL_TYPE::SMOOTH) {
      for (auto iVar = 0u; iVar < nVar; iVar++)
        nodes->SetSolution_Old(iPoint,iVar,0.0);

      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal). Covers all nVar (not just
            nu_tilde) so the Dirichlet solution[1..nVar-1]=0 imposed above for the Stochastic
            Backscatter Model's Langevin components (when active) is enforced consistently. ---*/

      if (implicit) {
        for (auto iVar = 0u; iVar < nVar; iVar++)
          Jacobian.DeleteValsRowi(iPoint, iVar);
      }
      continue;
    }

    /*--- For rough walls, the boundary condition is given by
    * (\frac{\partial \nu}{\partial n})_wall = \frac{\nu}{0.03*k_s}
    * where \nu is the solution variable, $n$ is the wall normal direction
    * and k_s is the equivalent sand grain roughness specified. ---*/

    /*--- Compute dual-grid area and boundary normal ---*/
    su2double Normal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);

    su2double Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Get laminar_viscosity and density ---*/
    su2double sigma = 2.0/3.0;
    su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double nu_lam = laminar_viscosity / density;

    su2double denom = sigma * 0.03 * Roughness_Height;
    su2double coeff = (nu_lam + nodes->GetSolution(iPoint,0)) / denom;

    su2double Res_Wall = Area * coeff * nodes->GetSolution(iPoint,0);
    LinSysRes(iPoint, 0) -= Res_Wall;

    su2double Jac_Wall = Area * (coeff + nodes->GetSolution(iPoint,0) / denom);
    if (implicit) Jacobian.AddVal2Diag(iPoint, 0, -Jac_Wall);
  }
  END_SU2_OMP_FOR
}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics*, CNumerics*,
                             CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];
  CFluidModel* FluidModel = flowSolver->GetFluidModel();
  const su2double* Turb_Properties = config->GetInlet_TurbVal(config->GetMarker_All_TagBound(val_marker));
  const su2double Nu_Factor = Turb_Properties[0];
  const su2double* Scalar_Inlet = config->GetKind_Species_Model() != SPECIES_MODEL::NONE
                                      ? config->GetInlet_SpeciesVal(config->GetMarker_All_TagBound(val_marker))
                                      : nullptr;

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    su2double nuTilde = Inlet_TurbVars[val_marker][iVertex][0];
    const auto* V_inlet = flowSolver->GetCharacPrimVar(val_marker, iVertex);

    /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
    if (config->GetInlet_Profile_From_File()) {
      nuTilde *= config->GetDensity_Ref() / config->GetViscosity_Ref();
    } else {
      /*--- Fluid model evaluation of the inlet nu tilde. ---*/
      su2double Density_Inlet;
      if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
        Density_Inlet = V_inlet[prim_idx.Density()];
        FluidModel->SetTDState_Prho(V_inlet[prim_idx.Pressure()], Density_Inlet);
      } else {
        FluidModel->SetTDState_T(V_inlet[prim_idx.Temperature()], Scalar_Inlet);
        Density_Inlet = FluidModel->GetDensity();
      }
      const su2double Laminar_Viscosity_Inlet = FluidModel->GetLaminarViscosity();
      nuTilde = Nu_Factor * Laminar_Viscosity_Inlet / Density_Inlet;
      if (config->GetSAParsedOptions().bc) nuTilde *= 0.005;
    }
    ghostNodes->SetSolution(iVertex, 0, nuTilde);

    /*--- The ghost rows are shared by all markers: the Langevin variables of stochastic backscatter
          (if any) must be set explicitly, they enter at zero. ---*/
    for (auto iVar = 1u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, 0.0);

    SetGhostPrimitives(iVertex, V_inlet);

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics*, CNumerics*,
                              CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Neumann: the turbulent variable is copied from the interior before computing the flux. ---*/
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nodes->GetSolution(iPoint, iVar));

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                     CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Neumann: the turbulent variable is copied from the interior before computing the flux. ---*/
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nodes->GetSolution(iPoint, iVar));

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                      CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    /*--- Prescribed turbulent state for an inflow. ---*/
    for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nu_tilde_Engine[iVar]);

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    SetGhostGeometry(geometry, val_marker, iVertex);
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config,  val_marker, true);
}

void CTurbSASolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config,  val_marker, false);
}

void CTurbSASolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container,
                               CNumerics*, CNumerics*,
                               CConfig *config, unsigned short val_marker, bool val_inlet_surface) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto GlobalIndex_donor = flowSolver->GetDonorGlobalIndex(val_marker, iVertex);
    const auto GlobalIndex = geometry->nodes->GetGlobalIndex(iPoint);

    /*--- No flux at a vertex whose donor is the point itself. ---*/
    if (GlobalIndex == GlobalIndex_donor) {
      ghostSkip[iVertex] = true;
      continue;
    }

    /*--- Normal vector for this vertex (negate for outward convention) ---*/

    su2double Normal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);

    su2double Area = GeometryToolbox::Norm(nDim, Normal);

    su2double UnitNormal[MAXNDIM] = {0.0};
    for (auto iDim = 0u; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    const auto* V_domain = flowSolver->GetNodes()->GetPrimitive(iPoint);

    /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/

    su2double Vn = GeometryToolbox::DotProduct(nDim, &V_domain[1], UnitNormal);

    bool ReverseFlow = false;
    if ((val_inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
    if ((!val_inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }

    /*--- No flux at all if there is a reverse flow, Euler b.c. for the direct problem. ---*/

    if (ReverseFlow) {
      ghostSkip[iVertex] = true;
      continue;
    }

    SetGhostPrimitives(iVertex, flowSolver->GetCharacPrimVar(val_marker, iVertex));

    /*--- Inflow analysis (interior extrapolation, a Neumann BC) or outflow analysis
     * (prescribed for an inflow), depending on which side val_marker is. ReverseFlow is always
     * false here (the other case returned above), so this reduces to val_inlet_surface. ---*/

    if (val_inlet_surface) {
      for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nodes->GetSolution(iPoint, iVar));
    } else {
      for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(iVertex, iVar, nu_tilde_ActDisk[iVar]);
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) ghostNormal(iVertex, iDim) = Normal[iDim];

    ghostSkip[iVertex] = false;
  }
  END_SU2_OMP_FOR

  /*--- The diffusive term at this boundary causes serious convergence problems. ---*/
  BoundaryFlux(geometry, solver_container, config,
               ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb()), val_marker);
}

void CTurbSASolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                         CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- The span loop below reaches a vertex of val_marker through GetOldVertex, which need not
   * cover every one of them, so every vertex starts skipped and only the ones actually filled
   * are cleared. ---*/
  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) ghostSkip[iVertex] = true;
  END_SU2_OMP_FOR

  for (auto iSpan = 0u; iSpan < nSpanWiseSections; iSpan++){
    su2double extAverageNu[MAXNVAR] = {0.0};
    extAverageNu[0] = flowSolver->GetMixingState(val_marker, iSpan, 5);

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(oldVertex, iVar, extAverageNu[iVar]);

      SetGhostPrimitives(oldVertex, flowSolver->GetCharacPrimVar(val_marker, oldVertex));

      SetGhostGeometry(geometry, val_marker, oldVertex);

      SetGhostDiffusionState(geometry, oldVertex, iPoint, Point_Normal);
    }
    END_SU2_OMP_FOR
  }

  BoundaryFlux(geometry, solver_container, config, ScalarFluxOptions::BoundaryFull(*config), val_marker);
}

void CTurbSASolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                   CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  auto* flowSolver = solver_container[FLOW_SOL];
  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  CFluidModel *FluidModel = flowSolver->GetFluidModel();

  su2double Factor_nu_Inf = config->GetNuFactor_FreeStream();

  /*--- The span loop below reaches a vertex of val_marker through GetOldVertex, which need not
   * cover every one of them, so every vertex starts skipped and only the ones actually filled
   * are cleared. ---*/
  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) ghostSkip[iVertex] = true;
  END_SU2_OMP_FOR

  /*--- Loop over all the spans on this boundary marker ---*/
  for (auto iSpan = 0; iSpan < nSpanWiseSections; iSpan++) {

    su2double rho       = flowSolver->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = flowSolver->GetAveragePressure(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double nu_tilde[MAXNVAR] = {0.0};
    nu_tilde[0] = Factor_nu_Inf*muLam/rho;

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      for (auto iVar = 0u; iVar < nVar; iVar++) ghostNodes->SetSolution(oldVertex, iVar, nu_tilde[iVar]);

      SetGhostPrimitives(oldVertex, flowSolver->GetCharacPrimVar(val_marker, oldVertex));

      SetGhostGeometry(geometry, val_marker, oldVertex);

      SetGhostDiffusionState(geometry, oldVertex, iPoint, Point_Normal);
    }
    END_SU2_OMP_FOR
  }

  BoundaryFlux(geometry, solver_container, config, ScalarFluxOptions::BoundaryFull(*config), val_marker);
}

void CTurbSASolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics*,
                                       CNumerics*, CConfig *config) {
  SU2_ZONE_SCOPED

  if (solver_container[FLOW_SOL] == nullptr) return;

  const auto optConv = ScalarFluxOptions::BoundaryConvective(*config, config->GetBounded_Turb());
  const auto optVisc = ScalarFluxOptions::BoundaryDiffusive(*config, true);

  /*--- SA's diffusion coefficients read no auxiliary ghost field. ---*/
  const auto fillGhostExtras = [](unsigned long, unsigned long) {};

  DispatchScheme<CScalarFlux_SA, 1, 4>(config, [&](auto tag) {
    FluidInterfaceFluxResidual<typename decltype(tag)::type>(geometry, solver_container, config, optConv, optVisc,
                                                             fillGhostExtras);
  });
}

void CTurbSASolver::SetTurbVars_WF(CGeometry *geometry, CSolver **solver_container,
                                  const CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- We use a very high max nr of iterations, but we only need this the first couple of iterations ---*/
  const unsigned short max_iter = config->GetwallModel_MaxIter();

  /* --- tolerance has LARGE impact on convergence, do not increase this value! --- */
  const su2double tol = 1e-12;


  /*--- Typical constants from boundary layer theory ---*/

  const su2double cv1_3 = 7.1*7.1*7.1;

  CVariable* flow_nodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint_Neighbor)) {

      su2double Y_Plus = solver_container[FLOW_SOL]->GetYPlus(val_marker, iVertex);

      /*--- Do not use wall model at the ipoint when y+ < "limit" ---*/

      if (Y_Plus < config->GetwallModel_MinYPlus()) continue;

      su2double Lam_Visc_Normal = flow_nodes->GetLaminarViscosity(iPoint_Neighbor);
      su2double Density_Normal = flow_nodes->GetDensity(iPoint_Neighbor);
      su2double Kin_Visc_Normal = Lam_Visc_Normal/Density_Normal;

      su2double Eddy_Visc = solver_container[FLOW_SOL]->GetEddyViscWall(val_marker, iVertex);

      /*--- Solve for the new value of nu_tilde given the eddy viscosity and using a Newton method.
       * f(nu_til) = nu_til^4 - nu_t*(nu_til^3 + cv1^3*nu^3) has a single positive root, which is larger than nu_t,
       * and f is increasing and convex for nu_til > nu_t, so the iteration converges from any start >= nu_t.
       * Start from the value imposed at the previous iteration, not from the wall node: with wall functions the
       * wall value is free and tends to zero (f' = 0) or below, e.g. after a restart without wall functions. ---*/

      const su2double nu_til_start = max(nodes->GetSolution(iPoint_Neighbor,0), Eddy_Visc/Density_Normal);
      su2double nu_til = 0.0;
      su2double nu_til_old = nu_til_start;

      unsigned short counter = 0;
      su2double diff = 1.0;
      su2double relax = config->GetwallModel_RelFac();
      while (diff > tol) {
        // note the error in Nichols and Nelson
        su2double func = pow(nu_til_old,4) - (Eddy_Visc/Density_Normal)*(pow(nu_til_old,3) + pow(Kin_Visc_Normal,3)*cv1_3);
        su2double func_prim = 4.0 * pow(nu_til_old,3) - 3.0*(Eddy_Visc/Density_Normal)*pow(nu_til_old,2);

        // damped Newton method
        nu_til = nu_til_old - relax*(func/func_prim);

        diff = fabs(nu_til-nu_til_old);
        nu_til_old = nu_til;

        // sometimes we get negative values when the solution has not converged yet, we just reset the nu_tilde in that case.
        if (nu_til_old<tol) {
          relax /= 2.0;
          nu_til_old = nu_til_start/relax;
        }

        counter++;
        if (counter > max_iter) break;
      }

      su2double nuTil[MAXNVAR] = {0.0};
      nuTil[0] = nu_til;
      nodes->SetSolution_Old(iPoint_Neighbor, nuTil);
      LinSysRes.SetBlock_Zero(iPoint_Neighbor);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal). Covers all nVar (not just
            nu_tilde) so the Dirichlet solution[1..nVar-1]=0 imposed above for the Stochastic
            Backscatter Model's Langevin components (when active) is enforced consistently. ---*/

      if (implicit) {
        for (auto iVar = 0u; iVar < nVar; iVar++)
          Jacobian.DeleteValsRowi(iPoint_Neighbor, iVar);
      }
    }
  }
}

void CTurbSASolver::SetDES_LengthScale(CSolver **solver, CGeometry *geometry, CConfig *config){
  SU2_ZONE_SCOPED

  const auto kindHybridRANSLES = config->GetKind_HybridRANSLES();

  const su2double constDES = config->GetConst_DES();

  const su2double fw_star = 0.424, cv1_3 = pow(7.1, 3), k2 = pow(0.41, 2);
  const su2double cb1   = 0.1355, ct3 = 1.2, ct4 = 0.5;
  const su2double sigma = 2./3., cb2 = 0.622, f_max = 1.0, f_min = 0.1, a1 = 0.15, a2 = 0.3;

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver[FLOW_SOL]->GetNodes());

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++){

    const auto coord_i       = geometry->nodes->GetCoord(iPoint);
    const auto nNeigh        = geometry->nodes->GetnPoint(iPoint);
    const auto wallDistance  = geometry->nodes->GetWall_Distance(iPoint);
    const auto velocityGrad  = flowNodes->GetVelocityGradient(iPoint);
    const auto vorticity     = flowNodes->GetVorticity(iPoint);
    const auto density       = flowNodes->GetDensity(iPoint);
    const auto laminarViscosity = flowNodes->GetLaminarViscosity(iPoint);
    const auto eddyViscosity    = nodes->GetmuT(iPoint);
    const su2double kinematicViscosity     = laminarViscosity/density;
    const su2double kinematicViscosityTurb = eddyViscosity/density;

    su2double uijuij = 0.0;
    for(auto iDim = 0u; iDim < nDim; iDim++){
      for(auto jDim = 0u; jDim < nDim; jDim++){
        uijuij += pow(velocityGrad[iDim][jDim], 2);
      }
    }
    uijuij = sqrt(fabs(uijuij));
    uijuij = max(uijuij,1e-10);

    /*--- Low Reynolds number correction term ---*/

    const su2double nu_hat = nodes->GetSolution(iPoint,0);
    const su2double Ji   = nu_hat/kinematicViscosity;
    const su2double Ji_2 = Ji * Ji;
    const su2double Ji_3 = Ji*Ji*Ji;
    const su2double fv1  = Ji_3/(Ji_3+cv1_3);
    const su2double fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    const su2double ft2 = ct3*exp(-ct4*Ji_2);
    const su2double cw1 = cb1/k2+(1.0+cb2)/sigma;

    su2double psi_2 = (1.0 - (cb1/(cw1*k2*fw_star))*(ft2 + (1.0 - ft2)*fv2))/(fv1 * max(1.0e-10,1.0-ft2));
    psi_2 = min(100.0,psi_2);

    su2double lengthScale = 0.0, lesSensor = 0.0, desFilterWidth = 0.0;

    const su2double LES_FilterWidth = config->GetLES_FilterWidth();

    switch(kindHybridRANSLES){
      case SA_DES: {
        /*--- Original Detached Eddy Simulation (DES97)
        Spalart
        1997
        ---*/

        desFilterWidth = geometry->nodes->GetMaxLength(iPoint);
        if (LES_FilterWidth > 0.0){
          desFilterWidth = LES_FilterWidth;
        }
        const su2double distDES = constDES * desFilterWidth;
        lengthScale = min(distDES,wallDistance);
        lesSensor = (wallDistance<=distDES) ? 0.0 : 1.0;

        if (config->GetEnforceLES()) {
          lengthScale = distDES;
          lesSensor = 1.0;
        }

        break;
      }
      case SA_DDES: {
        /*--- A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
         Spalart et al.
         Theoretical and Computational Fluid Dynamics - 2006
         ---*/

        desFilterWidth = geometry->nodes->GetMaxLength(iPoint);
        if (LES_FilterWidth > 0.0){
          desFilterWidth = LES_FilterWidth;
        }

        const su2double r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        const su2double f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        const su2double distDES = constDES * desFilterWidth;
        /*--- Same as wallDistance-f_d*max(0,wallDistance-distDES), without the cancellation of that
              difference when wallDistance >> distDES (e.g. no viscous walls, see CGeometry). ---*/
        lengthScale = (wallDistance<=distDES) ? su2double(wallDistance) : su2double((1.0-f_d)*wallDistance + f_d*distDES);
        lesSensor = (wallDistance<=distDES) ? 0.0 : f_d;

        if (config->GetEnforceLES()) {
          lengthScale = distDES;
          lesSensor = 1.0;
        }

        break;
      }
      case SA_ZDES: {
        /*--- Recent improvements in the Zonal Detached Eddy Simulation (ZDES) formulation.
         Deck
         Theoretical and Computational Fluid Dynamics - 2012
         ---*/

        const su2double deltaDDES = geometry->nodes->GetMaxLength(iPoint);

        su2double delta[MAXNDIM] = {}, ratioOmega[MAXNDIM] = {};

        for (const auto jPoint : geometry->nodes->GetPoints(iPoint)) {
          const auto coord_j = geometry->nodes->GetCoord(jPoint);
          for (auto iDim = 0u; iDim < nDim; iDim++){
            const su2double deltaAux = abs(coord_j[iDim] - coord_i[iDim]);
            delta[iDim] = max(delta[iDim], deltaAux);
          }
        }

        const su2double omega = GeometryToolbox::Norm(3, vorticity);
        const bool zeroVorticity = (omega <= EPS);

        for (auto iDim = 0u; iDim < 3; iDim++){
          ratioOmega[iDim] = zeroVorticity ? su2double(0.0) : su2double(vorticity[iDim]/omega);
        }

        desFilterWidth = sqrt(pow(ratioOmega[0], 2)*delta[1]*delta[2] +
                              pow(ratioOmega[1], 2)*delta[0]*delta[2] +
                              pow(ratioOmega[2], 2)*delta[0]*delta[1]);

        const su2double r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        const su2double f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        /*--- The vorticity-based filter width is undefined without vorticity. ---*/
        if (f_d < 0.99 || zeroVorticity){
          desFilterWidth = deltaDDES;
        }

        if (LES_FilterWidth > 0.0){
          desFilterWidth = LES_FilterWidth;
        }
        const su2double distDES = constDES * desFilterWidth;
        /*--- Same as wallDistance-f_d*max(0,wallDistance-distDES), without the cancellation of that
              difference when wallDistance >> distDES (e.g. no viscous walls, see CGeometry). ---*/
        lengthScale = (wallDistance<=distDES) ? su2double(wallDistance) : su2double((1.0-f_d)*wallDistance + f_d*distDES);
        lesSensor = (wallDistance<=distDES) ? 0.0 : f_d;

        if (config->GetEnforceLES()) {
          lengthScale = distDES;
          lesSensor = 1.0;
        }

        break;
      }
      case SA_EDDES: {
        /*--- An Enhanced Version of DES with Rapid Transition from RANS to LES in Separated Flows.
         Shur et al.
         Flow Turbulence Combust - 2015
         ---*/

        su2double vortexTiltingMeasure = nodes->GetVortex_Tilting(iPoint);

        const su2double omega = GeometryToolbox::Norm(3, vorticity);
        const bool zeroVorticity = (omega <= EPS);

        su2double ratioOmega[MAXNDIM] = {};

        for (auto iDim = 0; iDim < 3; iDim++){
          ratioOmega[iDim] = zeroVorticity ? su2double(0.0) : su2double(vorticity[iDim]/omega);
        }

        const su2double deltaDDES = geometry->nodes->GetMaxLength(iPoint);

        su2double ln_max = 0.0;
        for (const auto jPoint : geometry->nodes->GetPoints(iPoint)) {
          const auto coord_j = geometry->nodes->GetCoord(jPoint);
          su2double delta[MAXNDIM] = {};
          for (auto iDim = 0u; iDim < nDim; iDim++){
            delta[iDim] = fabs(coord_j[iDim] - coord_i[iDim]);
          }
          su2double ln[3];
          ln[0] = delta[1]*ratioOmega[2] - delta[2]*ratioOmega[1];
          ln[1] = delta[2]*ratioOmega[0] - delta[0]*ratioOmega[2];
          ln[2] = delta[0]*ratioOmega[1] - delta[1]*ratioOmega[0];
          const su2double aux_ln = sqrt(ln[0]*ln[0] + ln[1]*ln[1] + ln[2]*ln[2]);
          ln_max = max(ln_max, aux_ln);
          vortexTiltingMeasure += nodes->GetVortex_Tilting(jPoint);
        }
        vortexTiltingMeasure /= (nNeigh + 1);

        const su2double f_kh = max(f_min,
                                   min(f_max,
                                       f_min + ((f_max - f_min)/(a2 - a1)) * (vortexTiltingMeasure - a1)));

        const su2double r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
        const su2double f_d = 1.0-tanh(pow(8.0*r_d,3.0));

        desFilterWidth = (ln_max/sqrt(3.0)) * f_kh;
        /*--- The shear-layer-adapted filter width is undefined without vorticity. ---*/
        if (f_d < 0.999 || zeroVorticity){
          desFilterWidth = deltaDDES;
        }

        if (LES_FilterWidth > 0.0){
          desFilterWidth = LES_FilterWidth;
        }
        const su2double distDES = constDES * desFilterWidth;
        /*--- Same as wallDistance-f_d*max(0,wallDistance-distDES), without the cancellation of that
              difference when wallDistance >> distDES (e.g. no viscous walls, see CGeometry). ---*/
        lengthScale = (wallDistance<=distDES) ? su2double(wallDistance) : su2double((1.0-f_d)*wallDistance + f_d*distDES);
        lesSensor = (wallDistance<=distDES) ? 0.0 : f_d;

        if (config->GetEnforceLES()) {
          lengthScale = distDES;
          lesSensor = 1.0;
        }

        break;
      }
    }

    nodes->SetDES_LengthScale(iPoint, lengthScale);
    nodes->SetLES_Mode(iPoint, lesSensor);
    nodes->SetDES_FilterWidth(iPoint, desFilterWidth);

  }
  END_SU2_OMP_FOR
}

void CTurbSASolver::SetBackscatterInBox(CConfig *config, CGeometry *geometry) {
  SU2_ZONE_SCOPED

  auto sbsBoxBounds = config->GetSBSParam().StochBackscatterBoxBounds;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    const auto coord = geometry->nodes->GetCoord(iPoint);
    bool outOfBoxX = (coord[0]<sbsBoxBounds[0] || coord[0]>sbsBoxBounds[1]);
    bool outOfBoxY = (coord[1]<sbsBoxBounds[2] || coord[1]>sbsBoxBounds[3]);
    bool outOfBoxZ = (coord[2]<sbsBoxBounds[4] || coord[2]>sbsBoxBounds[5]);
    bool outOfBox  = (outOfBoxX || outOfBoxY || outOfBoxZ);
    su2double sbsInBox = outOfBox ? 0.0 : 1.0;
    nodes->SetSBSInBox(iPoint, sbsInBox);
  }
  END_SU2_OMP_FOR

}

void CTurbSASolver::SetLangevinSourceTerms(CConfig *config, CGeometry* geometry) {
  SU2_ZONE_SCOPED

  const su2double threshold = config->GetSBSParam().stochFdThreshold;
  const su2double dummySource = 1e3;
  unsigned long timeIter = config->GetTimeIter();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++){
    unsigned long iPointGlobal = geometry->nodes->GetGlobalIndex(iPoint);
    for (unsigned short iDim = 0; iDim < nDim; iDim++){
      su2double lesSensor = nodes->GetLES_Mode(iPoint) * nodes->GetSBSInBox(iPoint);
      if (lesSensor>threshold) {
        su2double rnd = RandomToolbox::GetNormal(iPointGlobal, iDim, timeIter);
        nodes->SetLangevinSourceTermsOld(iPoint, iDim, rnd);
        nodes->SetLangevinSourceTerms(iPoint, iDim, rnd);
      } else {
        nodes->SetLangevinSourceTermsOld(iPoint, iDim, dummySource);
        nodes->SetLangevinSourceTerms(iPoint, iDim, 0.0);
      }
    }
  }
  END_SU2_OMP_FOR

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto kindBC = config->GetMarker_All_KindBC(iMarker);
      if (kindBC != SEND_RECEIVE && kindBC != PERIODIC_BOUNDARY) {
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          nodes->SetLangevinSourceTermsOld(iPoint, iDim, dummySource);
          nodes->SetLangevinSourceTerms(iPoint, iDim, 0.0);
        }
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- The two copies of a periodic point have different global indices, hence would draw different
        random numbers: copy the values (and the active/inactive flag, encoded in the "old" value)
        from the master to the passive periodic face. ---*/

  if (config->GetnMarker_Periodic() > 0) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSBSPeriodicBuf(iPoint, iDim, nodes->GetLangevinSourceTermsOld(iPoint, iDim));
        nodes->SetSBSPeriodicBuf(iPoint, 3+iDim, nodes->GetLangevinSourceTerms(iPoint, iDim));
      }
    }
    END_SU2_OMP_FOR

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_COPY);

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        nodes->SetLangevinSourceTermsOld(iPoint, iDim, nodes->GetSBSPeriodicBuf(iPoint, iDim));
        nodes->SetLangevinSourceTerms(iPoint, iDim, nodes->GetSBSPeriodicBuf(iPoint, 3+iDim));
      }
    }
    END_SU2_OMP_FOR
  }
}

void CTurbSASolver::ComputeOU_Process(CSolver **solver, CConfig *config, CGeometry *geometry) {
  SU2_ZONE_SCOPED

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver[FLOW_SOL]->GetNodes());
  su2double timeStep = config->GetDelta_UnstTimeND();

  /*--- Same time scale as the Langevin equations (CSourceBase_TurbSA::ResidualStochEquations): LES time
        scale Ctau*Delta^2/nu_t blended by the LES sensor with the (bounded) RANS one d_w^2/nu_t, and the
        time-averaged eddy viscosity if SBS_USE_MEAN_TURB, consistently with the SST counterpart. ---*/
  const bool useMeanTurb = config->GetSBSParam().useMeanTurb;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double maxDelta = nodes->GetDES_FilterWidth(iPoint);
    su2double wallDist = geometry->nodes->GetWall_Distance(iPoint);
    su2double nuT = useMeanTurb ? nodes->GetMeanEddyViscosity(iPoint)
                                : nodes->GetmuT(iPoint) / flowNodes->GetDensity(iPoint);
    nuT = max(nuT, 1e-10);
    su2double timeLES = config->GetSBSParam().SBS_Ctau * maxDelta*maxDelta / nuT;
    su2double timeRANS = min(wallDist*wallDist / nuT, 10.0*timeStep);
    su2double lesMode = nodes->GetLES_Mode(iPoint);
    su2double turbTime = timeRANS * (1.0-lesMode) + timeLES * lesMode;
    su2double timeRatio = timeStep/turbTime;
    su2double term1 = exp(-timeRatio);
    su2double term2 = (timeRatio < 1e-6) ? sqrt(2.0*timeRatio) : sqrt(1.0-exp(-2.0*timeRatio));
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      su2double stochSourceOld = nodes->GetOU_Process(iPoint, iDim);
      su2double langevinSource = nodes->GetLangevinSourceTerms(iPoint, iDim);
      su2double stochSourceNew = stochSourceOld*term1 + term2*langevinSource;
      nodes->SetOU_Process(iPoint, iDim, stochSourceNew);
    }
  }
  END_SU2_OMP_FOR

  /*--- Keep both copies of a periodic point identical (the time scale depends on local quantities,
        e.g. the filter width, which can differ slightly on the two sides). ---*/

  if (config->GetnMarker_Periodic() > 0) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        nodes->SetSBSPeriodicBuf(iPoint, iDim, nodes->GetOU_Process(iPoint, iDim));
    END_SU2_OMP_FOR

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_COPY);

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        nodes->SetOU_Process(iPoint, iDim, nodes->GetSBSPeriodicBuf(iPoint, iDim));
    END_SU2_OMP_FOR
  }

  InitiateComms(geometry, config, MPI_QUANTITIES::OU_PROCESS);
  CompleteComms(geometry, config, MPI_QUANTITIES::OU_PROCESS);
}

void CTurbSASolver::SmoothLangevinSourceTerms(CConfig* config, CGeometry* geometry) {
  SU2_ZONE_SCOPED

  static std::array<su2double, 3> globalResNorm;
  static std::array<unsigned long, 3> global_nPointLES;
  static std::array<su2double, 6> globalChecks;
  /*--- Scalars (one per component) driving BiCGSTAB control flow (loop conditions/break). These
        must be "static" (as in the original Jacobi implementation) rather than plain locals: the
        whole point-loop and Allreduce machinery below runs master-thread-only per MPI rank (BEGIN/
        END_SU2_OMP_SAFE_GLOBAL_ACCESS), but the loop conditions and breakdown checks are plain code
        executed redundantly by every OpenMP thread sharing this parallel region, so they must all
        observe the exact same, already-synchronized value that only the master thread computed. ---*/
  static std::array<su2double, 3> globalRho, globalRhoPrev, globalAlpha, globalOmega, globalR0V, globalTS, globalTT;

  const su2double cDelta = config->GetSBSParam().SBS_Cdelta;
  const unsigned short maxIter = config->GetSBSParam().SBS_maxIterSmooth;
  const su2double tol = -5.0;
  const su2double sourceLim = 5.0;
  const su2double eps_breakdown = 1e-50;
  /*--- The non-orthogonal correction only needs to be refreshed periodically (its role is to
        accelerate/correct convergence to the true gradient-based solution, not to change the fixed
        point), so its extra point loop and MPI round-trip are skipped on iterations in between. It
        also defines the RHS of the linear system solved below, which must stay fixed while BiCGSTAB
        builds its Krylov subspace, so a refresh restarts (warm-started from the current solution)
        the Krylov iteration. ---*/
  const unsigned short gradRefreshInterval = 5;
  unsigned long timeIter = config->GetTimeIter();
  unsigned long restartIter = config->GetRestart_Iter();

  /*--- Assemble system matrix: the orthogonal (implicit) coefficient a_ij and the diagonal
        diag_i = 1 + sum(a_ij), both purely geometric (independent of iDim), plus the
        non-orthogonal correction vector betaVec used to reconstruct the full gradient-based
        diffusive flux across each face (deferred correction, folded into the RHS, see below).
        All edge-based sums (diagonal, matrix-vector products, gradient, non-orthogonal RHS) are
        partial at periodic points, since each copy of the point only sees the edges on its own
        side: they are completed with the contribution of the matching point (PERIODIC_SBS_SUM)
        before being used, so that both copies hold the same (complete) values. */

  if (timeIter == restartIter) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      su2double maxDelta = nodes->GetDES_FilterWidth(iPoint);
      su2double b2 = cDelta * maxDelta * maxDelta;
      su2double volume_iPoint = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      auto coord_i = geometry->nodes->GetCoord(iPoint);
      const auto nNeighbors = geometry->nodes->GetnPoint(iPoint);

      /*--- smoothMatrix/smoothBetaVec (see CTurbSAVariable) are fixed-size, MAXNNEIGHBORS wide;
            silently exceeding that bound here would overrun those buffers. ---*/
      if (nNeighbors > CTurbSAVariable::MAXNNEIGHBORS) {
        SU2_MPI::Error("Point " + std::to_string(iPoint) + " has " + std::to_string(nNeighbors) +
                       " point-to-point neighbors, exceeding CTurbSAVariable::MAXNNEIGHBORS (" +
                       std::to_string(CTurbSAVariable::MAXNNEIGHBORS) + "). Increase MAXNNEIGHBORS "
                       "in CTurbSAVariable.hpp and recompile.", CURRENT_FUNCTION);
      }

      su2double sumA = 0.0;
      su2double betaTensor[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // xx, yy, zz, xy, yz, xz
      for (unsigned short iNode = 0; iNode < nNeighbors; iNode++) {
        auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
        auto coord_j = geometry->nodes->GetCoord(jPoint);
        auto iEdge = geometry->nodes->GetEdge(iPoint, iNode);
        auto* normal = geometry->edges->GetNormal(iEdge);
        su2double area = GeometryToolbox::Norm(nDim, normal);
        su2double dx_ij_vec[3];
        for (unsigned short index = 0; index < nDim; index++)
          dx_ij_vec[index] = coord_j[index] - coord_i[index];
        su2double distance = GeometryToolbox::Norm(nDim, dx_ij_vec);
        su2double dist_ij_2 = max(distance*distance, 1e-10);
        su2double dot_nd = 0.0;
        for (unsigned short index = 0; index < nDim; index++)
          dot_nd += normal[index] * dx_ij_vec[index];
        su2double sign = (geometry->edges->GetNode(iEdge, 0) == iPoint) ? 1.0 : -1.0;
        su2double d_normal = sign * dot_nd / max(area*distance, 1e-10);
        su2double a_ij = area/volume_iPoint * fabs(d_normal) * b2/max(distance, 1e-10);
        nodes->SetSmoothingMatrixCoeff(iPoint, iNode, a_ij);
        sumA += a_ij;

        /*--- Directional decomposition of the operator, 0.5*sum_j(a_ij * e_ij e_ij^T) with e_ij the
              unit edge vector: on a Cartesian grid its eigenvalues are exactly the coefficients
              beta_x, beta_y, beta_z of the equivalent lattice operator assumed by the Bessel
              scaling (diag = 1 + 2*(beta_x+beta_y+beta_z)), consistently with a_ij by construction. ---*/

        const su2double halfA = 0.5 * a_ij / dist_ij_2;
        betaTensor[0] += halfA * dx_ij_vec[0] * dx_ij_vec[0];
        betaTensor[1] += halfA * dx_ij_vec[1] * dx_ij_vec[1];
        betaTensor[2] += halfA * dx_ij_vec[2] * dx_ij_vec[2];
        betaTensor[3] += halfA * dx_ij_vec[0] * dx_ij_vec[1];
        betaTensor[4] += halfA * dx_ij_vec[1] * dx_ij_vec[2];
        betaTensor[5] += halfA * dx_ij_vec[0] * dx_ij_vec[2];

        /*--- betaVec = (b2/V_i) * sign * (normal - (dot_nd/dist_ij_2)*edge_vector), the coefficient
              such that mean_grad_face . betaVec gives the non-orthogonal (tangential) part of the
              gradient-based diffusive flux, consistent with CAvgGrad_Base::CorrectGradient. ---*/

        for (unsigned short index = 0; index < nDim; index++) {
          su2double betaVec = (b2/volume_iPoint) * sign * (normal[index] - (dot_nd/dist_ij_2) * dx_ij_vec[index]);
          nodes->SetSmoothingBetaVec(iPoint, iNode, index, betaVec);
        }
      }
      nodes->SetSBSPeriodicBuf(iPoint, 0, sumA);
      for (unsigned short iComp = 0; iComp < 6; iComp++)
        nodes->SetSBSPeriodicBuf(iPoint, 1+iComp, betaTensor[iComp]);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_SUM);

    /*--- Diagonal coefficient and, if requested, the Bessel integral used to rescale the smoothed
          field to unit variance. The latter is purely geometric too, so it is computed once here
          (not once per component and per time step). ---*/

    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      nodes->SetSmoothingDiag(iPoint, 1.0 + nodes->GetSBSPeriodicBuf(iPoint, 0));

      if (!config->GetSBSParam().besselScaleFactor) continue;

      /*--- Eigenvalues (beta_I, beta_J, beta_K) of the symmetric 3x3 tensor betaTensor. ---*/

      su2double M[3][3];
      M[0][0] = nodes->GetSBSPeriodicBuf(iPoint, 1);
      M[1][1] = nodes->GetSBSPeriodicBuf(iPoint, 2);
      M[2][2] = nodes->GetSBSPeriodicBuf(iPoint, 3);
      M[0][1] = M[1][0] = nodes->GetSBSPeriodicBuf(iPoint, 4);
      M[1][2] = M[2][1] = nodes->GetSBSPeriodicBuf(iPoint, 5);
      M[0][2] = M[2][0] = nodes->GetSBSPeriodicBuf(iPoint, 6);

      su2double lambda[3] = {0.0};
      const su2double p1 = M[0][1]*M[0][1] + M[1][2]*M[1][2] + M[0][2]*M[0][2];
      if (p1 < 1e-20) {
        lambda[0] = M[0][0];
        lambda[1] = M[1][1];
        lambda[2] = M[2][2];
      } else {
        const su2double trace = (M[0][0] + M[1][1] + M[2][2]) / 3.0;
        const su2double p2 = pow(M[0][0]-trace, 2) + pow(M[1][1]-trace, 2) + pow(M[2][2]-trace, 2) + 2.0 * p1;
        const su2double p = sqrt(p2 / 6.0);
        su2double B[3][3];
        for (unsigned short ind1 = 0; ind1 < 3; ind1++)
          for (unsigned short ind2 = 0; ind2 < 3; ind2++)
            B[ind1][ind2] = (M[ind1][ind2] - ((ind1 == ind2) ? trace : 0.0)) / p;
        const su2double detB =
            B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[2][1]) -
            B[0][1]*(B[1][0]*B[2][2] - B[1][2]*B[2][0]) +
            B[0][2]*(B[1][0]*B[2][1] - B[1][1]*B[2][0]);
        const su2double r = max(min(0.5 * detB, 1.0), -1.0);
        const su2double phi = acos(r) / 3.0;
        lambda[0] = trace + 2.0*p*cos(phi);
        lambda[1] = trace + 2.0*p*cos(phi + 2.0*PI_NUMBER/3.0);
        lambda[2] = trace + 2.0*p*cos(phi + 4.0*PI_NUMBER/3.0);
      }
      for (unsigned short iComp = 0; iComp < 3; iComp++) lambda[iComp] = max(lambda[iComp], 0.0);

      nodes->SetBesselIntegral(iPoint, RandomToolbox::GetBesselIntegral(lambda[0], lambda[1], lambda[2]));
    }
    END_SU2_OMP_FOR
  }

  /*--- Points actually solved for (active, see SetLangevinSourceTerms), and points whose Green-Gauss
        gradient or zeroed phat/shat can be read by an active point (on this or another rank, or across
        a periodic boundary), per component, so that the point loops of the solver only visit them
        instead of the whole partition. The rank-private BiCGSTAB work vectors are sized to the active
        points. All of these are filled and used by the master thread only, like the point loops
        themselves (the other threads keep them empty). ---*/

  using WorkVecs = std::array<std::vector<su2double>, 3>;
  std::array<std::vector<unsigned long>, 3> activePoints, gradPoints;
  std::array<std::vector<bool>, 3> isActive;
  std::vector<bool> isPeriodic;
  WorkVecs bRhs, rVec, rhat0Vec, pVec, vVec, sVec, tVec;

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    isPeriodic.assign(nPointDomain, false);
    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) continue;
      for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (iPoint < nPointDomain) isPeriodic[iPoint] = true;
      }
    }

    unsigned long local_nPointLES[3] = {0, 0, 0};
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      isActive[iDim].assign(nPointDomain, false);
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
        isActive[iDim][iPoint] = true;
        activePoints[iDim].push_back(iPoint);
      }
      /*--- Halo neighbors are conservatively taken as possibly active (their flag is not exchanged),
            and periodic points always need their partial gradient for the periodic sum. ---*/
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        bool needed = isActive[iDim][iPoint] || isPeriodic[iPoint];
        for (unsigned short iNode = 0; !needed && iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
          const auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
          needed = (jPoint >= nPointDomain) || isActive[iDim][jPoint];
        }
        if (needed) gradPoints[iDim].push_back(iPoint);
      }
      for (auto* vec : {&bRhs, &rVec, &rhat0Vec, &pVec, &vVec, &sVec, &tVec}) (*vec)[iDim].assign(activePoints[iDim].size(), 0.0);
      local_nPointLES[iDim] = activePoints[iDim].size();
    }
    SU2_MPI::Allreduce(local_nPointLES, global_nPointLES.data(), nDim, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Matrix-free operator A(x)_i = diag_i*x_i - sum_j(a_ij*x_j), the only part of the discrete
        system that is a genuine matrix (the tangential/non-orthogonal term below is a fixed RHS
        contribution between gradient refreshes, not part of A), applied to the components in comps
        at once. GetInput(j, iDim) reads whichever field is currently being multiplied by A (the
        solution itself to form the initial residual, or the preconditioned BiCGSTAB directions
        phat/shat); the result is written into the local, rank-private (never used as a neighbor, so
        never halo-exchanged) work vectors. Slot iDim of the periodic buffer holds the partial sum of
        component iDim, so a single periodic exchange serves all the components. Runs master-thread-
        only per rank, exactly like the original Jacobi sweep. ---*/

  auto ComputeMatVec = [&](const std::vector<unsigned short>& comps, auto&& GetInput, WorkVecs& out) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iDim : comps) {
      for (const auto iPoint : activePoints[iDim]) {
        su2double sum = 0.0;
        for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
          auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
          su2double a_ij = nodes->GetSmoothingMatrixCoeff(iPoint, iNode);
          sum += a_ij * GetInput(jPoint, iDim);
        }
        nodes->SetSBSPeriodicBuf(iPoint, iDim, sum);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_SUM);

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iDim : comps) {
      for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
        const auto iPoint = activePoints[iDim][k];
        out[iDim][k] = nodes->GetSmoothingDiag(iPoint)*GetInput(iPoint, iDim) - nodes->GetSBSPeriodicBuf(iPoint, iDim);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  };

  /*--- Dot products between rank-private work vectors (a . b into out, for each job), restricted to
        the points actually being solved for (points outside the active LES region are frozen at 0 and
        never enter the linear system, mirroring the exclusion in the original Jacobi sweep), for the
        components in comps, all reduced over the MPI ranks in a single call. Writes into the static
        scalars "out" so every OpenMP thread sharing this parallel region observes the same, already-
        synchronized values once outside the master-only block. ---*/

  struct DotJob { const WorkVecs* a; const WorkVecs* b; std::array<su2double, 3>* out; };

  auto DotActive = [&](const std::vector<unsigned short>& comps, std::initializer_list<DotJob> jobs) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      su2double local[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, global[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      int nValues = 0;
      for (const auto& job : jobs) {
        for (const auto iDim : comps) {
          const auto& a = (*job.a)[iDim];
          const auto& b = (*job.b)[iDim];
          for (unsigned long k = 0; k < a.size(); k++) local[nValues] += a[k]*b[k];
          nValues++;
        }
      }
      SU2_MPI::Allreduce(local, global, nValues, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      nValues = 0;
      for (const auto& job : jobs)
        for (const auto iDim : comps) (*job.out)[iDim] = global[nValues++];
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  };

  /*--- Refresh the RHS of component iDim at the start (or restart) of one of its blocks: recompute the
        Green-Gauss gradient used for the non-orthogonal RHS correction (one component at a time, since
        the gradient storage holds a single component). The halo values of the solution must be up to
        date. ---*/

  auto RefreshRhs = [&](unsigned short iDim) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iPoint : gradPoints[iDim]) {
      su2double source_i = nodes->GetLangevinSourceTerms(iPoint, iDim);
      su2double grad_i[3] = {0.0, 0.0, 0.0};
      for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
        auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
        auto iEdge = geometry->nodes->GetEdge(iPoint, iNode);
        auto* normal = geometry->edges->GetNormal(iEdge);
        su2double sign = (geometry->edges->GetNode(iEdge, 0) == iPoint) ? 1.0 : -1.0;
        su2double source_j = nodes->GetLangevinSourceTerms(jPoint, iDim);
        su2double phi_face = 0.5*(source_i + source_j);
        for (unsigned short index = 0; index < nDim; index++)
          grad_i[index] += sign * phi_face * normal[index];
      }
      for (unsigned short index = 0; index < nDim; index++)
        nodes->SetSBSPeriodicBuf(iPoint, index, grad_i[index]);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- At periodic points the edge-based surface integral is completed with the edges of the
          matching point; the (equal and opposite) periodic boundary faces of the two copies,
          omitted on both sides, cancel out exactly. ---*/

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_SUM);

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iPoint : gradPoints[iDim]) {
      su2double volume_iPoint = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      for (unsigned short index = 0; index < nDim; index++)
        nodes->SetLangevinSourceGrad(iPoint, index, nodes->GetSBSPeriodicBuf(iPoint, index) / max(volume_iPoint, 1e-10));
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG_GRAD);
    CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG_GRAD);

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      /*--- Inactive points are frozen at 0 (homogeneous Dirichlet) and must enter the matrix-free
            operator as zeros: clear any phat/shat left over from a previous time step in which
            the point was active, at the inactive points that ComputeMatVec can read as neighbors. ---*/
      for (const auto iPoint : gradPoints[iDim]) {
        if (isActive[iDim][iPoint]) continue;
        nodes->SetSmoothPhat(iPoint, iDim, 0.0);
        nodes->SetSmoothShat(iPoint, iDim, 0.0);
      }

      for (const auto iPoint : activePoints[iDim]) {
        su2double tangential_i = 0.0;
        for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
          auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
          su2double tangential_ij = 0.0;
          for (unsigned short index = 0; index < nDim; index++) {
            su2double mean_grad = 0.5*(nodes->GetLangevinSourceGrad(iPoint, index) + nodes->GetLangevinSourceGrad(jPoint, index));
            tangential_ij += mean_grad * nodes->GetSmoothingBetaVec(iPoint, iNode, index);
          }
          tangential_i += tangential_ij;
        }
        nodes->SetSBSPeriodicBuf(iPoint, 0, tangential_i);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    SBSPeriodicComm(geometry, config, PERIODIC_SBS_SUM);

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
      const auto iPoint = activePoints[iDim][k];
      bRhs[iDim][k] = nodes->GetLangevinSourceTermsOld(iPoint, iDim) + nodes->GetSBSPeriodicBuf(iPoint, 0);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  };

  /*--- Solve, for each spatial dimension, diag_i*x_i - sum_j(a_ij*x_j) = source_i_old + tangential_i
        with a matrix-free, Jacobi(diagonal)-preconditioned BiCGSTAB, restarting the Krylov subspace
        (warm-started from the current solution) every time the tangential/gradient RHS term is
        refreshed, so the matrix and RHS are always fixed for the duration of a Krylov subspace.
        The components are independent systems, solved in lockstep: each one keeps its own BiCGSTAB
        scalars, restarts and stopping criteria (so it follows exactly the same iterations as if it
        were solved on its own), while the halo exchanges (which carry all the components anyway),
        the periodic exchanges and the global reductions of a step are shared by all the components
        still iterating. The control flow below is executed redundantly by every OpenMP thread, and
        only depends on the (static) globally reduced scalars, so all threads take the same path. ---*/

  std::array<unsigned short, 3> totalIter{}, printIter{}, blockIter{};
  std::array<bool, 3> converged{}, needRestart{};
  for (unsigned short iDim = 0; iDim < nDim; iDim++) needRestart[iDim] = true;

  /*--- Convergence history of each component, printed once all of them are done so that the
        output stays grouped per component (filled by the master thread of the master rank). ---*/
  std::array<std::vector<std::pair<unsigned short, su2double>>, 3> history;

  std::vector<unsigned short> comps;

  auto IsDone = [&](unsigned short iDim) { return converged[iDim] || totalIter[iDim] >= maxIter; };

  /*--- Components whose BiCGSTAB breaks down (denominator ~0) leave the current step and restart
        their block, exactly as when solved on their own. ---*/
  auto DropBrokenDown = [&](const std::array<su2double, 3>& denom) {
    std::vector<unsigned short> kept;
    for (const auto iDim : comps) {
      if (fabs(denom[iDim]) < eps_breakdown) needRestart[iDim] = true;
      else kept.push_back(iDim);
    }
    comps = kept;
  };

  while (true) {

    /*--- Start (or restart) of a block for the components that need it: refresh the halo values of
          the solution (all components at once), the RHS, and rebuild the initial residual
          r0 = b - A*x (warm-started from the current solution estimate). ---*/

    comps.clear();
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      if (!IsDone(iDim) && needRestart[iDim]) comps.push_back(iDim);

    if (!comps.empty()) {

      InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
      CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);

      for (const auto iDim : comps) RefreshRhs(iDim);

      ComputeMatVec(comps, [&](unsigned long j, unsigned short iDim){ return nodes->GetLangevinSourceTerms(j, iDim); }, rVec);

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      for (const auto iDim : comps) {
        for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
          rVec[iDim][k] = bRhs[iDim][k] - rVec[iDim][k];
          rhat0Vec[iDim][k] = rVec[iDim][k];
          /*--- p_0 = v_0 = 0 (standard BiCGSTAB init): with rho_prev=alpha=omega=1 reset below, the
                first inner iteration's p update (p = r + beta*(p-omega*v)) then correctly reduces to
                p_1 = r_0 regardless of beta, since the (p_0 - omega*v_0) term vanishes. ---*/
          pVec[iDim][k] = 0.0;
          vVec[iDim][k] = 0.0;
        }
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      DotActive(comps, {{&rVec, &rVec, &globalResNorm}});

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      for (const auto iDim : comps) {
        globalResNorm[iDim] = (global_nPointLES[iDim]==0) ? su2double(0.0) : sqrt(globalResNorm[iDim] / global_nPointLES[iDim]);
        globalRho[iDim] = 1.0; globalAlpha[iDim] = 1.0; globalOmega[iDim] = 1.0;
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      for (const auto iDim : comps) {
        /*--- Account for the cost of this refresh/restart against the iteration budget: this
              guarantees the block/restart loop always makes progress towards maxIter, even in the
              degenerate case where every restarted Krylov subspace breaks down immediately (e.g. an
              already-converged residual, for which rho would be ~0 by construction). ---*/
        totalIter[iDim]++;
        if (log10(globalResNorm[iDim]) < tol) converged[iDim] = true;
        blockIter[iDim] = 0;
        needRestart[iDim] = false;
      }
    }

    /*--- One BiCGSTAB step for all the components still iterating. After the restarts above, a
          component that is not done never needs a restart here, so none left means all done. ---*/

    comps.clear();
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      if (!IsDone(iDim) && !needRestart[iDim]) comps.push_back(iDim);

    if (comps.empty()) break;

    SU2_OMP_SAFE_GLOBAL_ACCESS(for (const auto iDim : comps) globalRhoPrev[iDim] = globalRho[iDim];)
    DotActive(comps, {{&rhat0Vec, &rVec, &globalRho}});
    DropBrokenDown(globalRho);
    if (comps.empty()) continue;

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iDim : comps) {
      const su2double beta = (globalRho[iDim]/globalRhoPrev[iDim]) * (globalAlpha[iDim]/globalOmega[iDim]);
      for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
        const auto iPoint = activePoints[iDim][k];
        pVec[iDim][k] = rVec[iDim][k] + beta*(pVec[iDim][k] - globalOmega[iDim]*vVec[iDim][k]);
        nodes->SetSmoothPhat(iPoint, iDim, pVec[iDim][k] / nodes->GetSmoothingDiag(iPoint));
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    InitiateComms(geometry, config, MPI_QUANTITIES::SMOOTH_PHAT);
    CompleteComms(geometry, config, MPI_QUANTITIES::SMOOTH_PHAT);

    ComputeMatVec(comps, [&](unsigned long j, unsigned short iDim){ return nodes->GetSmoothPhat(j, iDim); }, vVec);

    DotActive(comps, {{&rhat0Vec, &vVec, &globalR0V}});
    DropBrokenDown(globalR0V);
    if (comps.empty()) continue;
    SU2_OMP_SAFE_GLOBAL_ACCESS(for (const auto iDim : comps) globalAlpha[iDim] = globalRho[iDim] / globalR0V[iDim];)

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (const auto iDim : comps) {
      for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
        const auto iPoint = activePoints[iDim][k];
        sVec[iDim][k] = rVec[iDim][k] - globalAlpha[iDim]*vVec[iDim][k];
        nodes->SetSmoothShat(iPoint, iDim, sVec[iDim][k] / nodes->GetSmoothingDiag(iPoint));
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    InitiateComms(geometry, config, MPI_QUANTITIES::SMOOTH_SHAT);
    CompleteComms(geometry, config, MPI_QUANTITIES::SMOOTH_SHAT);

    ComputeMatVec(comps, [&](unsigned long j, unsigned short iDim){ return nodes->GetSmoothShat(j, iDim); }, tVec);

    DotActive(comps, {{&tVec, &sVec, &globalTS}, {&tVec, &tVec, &globalTT}});
    DropBrokenDown(globalTT);
    if (comps.empty()) continue;
    SU2_OMP_SAFE_GLOBAL_ACCESS(for (const auto iDim : comps) globalOmega[iDim] = globalTS[iDim] / globalTT[iDim];)

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      su2double localResNorm[3] = {0.0, 0.0, 0.0}, globalResNormSq[3] = {0.0, 0.0, 0.0};
      for (size_t iComp = 0; iComp < comps.size(); iComp++) {
        const auto iDim = comps[iComp];
        for (unsigned long k = 0; k < activePoints[iDim].size(); k++) {
          const auto iPoint = activePoints[iDim][k];
          su2double source_i = nodes->GetLangevinSourceTerms(iPoint, iDim);
          su2double phat = nodes->GetSmoothPhat(iPoint, iDim);
          su2double shat = nodes->GetSmoothShat(iPoint, iDim);
          source_i += globalAlpha[iDim]*phat + globalOmega[iDim]*shat;
          nodes->SetLangevinSourceTerms(iPoint, iDim, source_i);
          rVec[iDim][k] = sVec[iDim][k] - globalOmega[iDim]*tVec[iDim][k];
          localResNorm[iComp] += pow(rVec[iDim][k], 2);
        }
      }
      SU2_MPI::Allreduce(localResNorm, globalResNormSq, static_cast<int>(comps.size()), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      for (size_t iComp = 0; iComp < comps.size(); iComp++) {
        const auto iDim = comps[iComp];
        globalResNorm[iDim] = (global_nPointLES[iDim]==0) ? su2double(0.0) : sqrt(globalResNormSq[iComp] / global_nPointLES[iDim]);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    for (const auto iDim : comps) {
      blockIter[iDim]++; totalIter[iDim]++; printIter[iDim]++;

      /*--- Record convergence at the first iteration, then every 10, and finally once more at
            whichever iteration ends the solve (converged or maxIter reached). printIter counts only
            real BiCGSTAB update steps (unlike totalIter, which also advances once per block restart
            to bound the iterations), so the first-iteration record is the first actual step. ---*/
      const bool isConverged = (log10(globalResNorm[iDim]) < tol || totalIter[iDim] == maxIter);
      const bool printConvergence = (printIter[iDim] == 1) || (printIter[iDim] % 10 == 0) || isConverged;

      if (printConvergence) {
        SU2_OMP_MASTER
        if (rank == MASTER_NODE) history[iDim].emplace_back(printIter[iDim], log10(globalResNorm[iDim]));
        END_SU2_OMP_MASTER
      }

      converged[iDim] = isConverged;
      if (blockIter[iDim] >= gradRefreshInterval) needRestart[iDim] = true;
    }
  }

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {

    SU2_OMP_MASTER
    if (rank == MASTER_NODE) {
      cout << "\nResidual of Laplacian smoothing along dimension " << iDim+1
           << "\n---------------------------------"
           << "\n   Iter       RMS Residual"
           << "\n---------------------------------" << endl;
      for (const auto& record : history[iDim]) {
        cout << "  "
             << std::setw(5) << record.first
             << "       "
             << std::setw(12) << std::fixed << std::setprecision(6) << record.second
             << endl;
      }
      cout << "---------------------------------" << endl;
    }
    END_SU2_OMP_MASTER

    {

        /*--- Scale source terms for variance preservation. ---*/

        su2double mean_check_old = 0.0, var_check_old = 0.0;
        su2double mean_check_new = 0.0, var_check_new = 0.0;
        su2double mean_check_notSmoothed = 0.0, var_check_notSmoothed = 0.0;

        SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          su2double source_notSmoothed = nodes->GetLangevinSourceTermsOld(iPoint, iDim);
          if (source_notSmoothed > 3.0*sourceLim) continue;
          su2double source = nodes->GetLangevinSourceTerms(iPoint, iDim);
          mean_check_old += source;
          var_check_old += pow(source, 2);
          mean_check_notSmoothed += source_notSmoothed;
          var_check_notSmoothed += pow(source_notSmoothed, 2);
        }
        END_SU2_OMP_FOR

        if (config->GetSBSParam().besselScaleFactor) {
          SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
          for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
            /*--- Computed once during the matrix assembly, see above. ---*/
            const su2double integral = nodes->GetBesselIntegral(iPoint);
            su2double scaleFactor = 1.0 / sqrt(max(integral, 1e-10));
            su2double source = nodes->GetLangevinSourceTerms(iPoint, iDim);
            source *= scaleFactor;
            if (source < -sourceLim || source > sourceLim) source = 0.0;
            mean_check_new += source;
            var_check_new += pow(source, 2);
            nodes->SetLangevinSourceTerms(iPoint, iDim, source);
          }
          END_SU2_OMP_FOR
        }

        SU2_OMP_SAFE_GLOBAL_ACCESS(globalChecks = {0, 0, 0, 0, 0, 0};)

        atomicAdd(mean_check_old, globalChecks[0]);
        atomicAdd(var_check_old, globalChecks[1]);
        atomicAdd(mean_check_notSmoothed, globalChecks[2]);
        atomicAdd(var_check_notSmoothed, globalChecks[3]);
        atomicAdd(mean_check_new, globalChecks[4]);
        atomicAdd(var_check_new, globalChecks[5]);

        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
          auto tmp = globalChecks;
          SU2_MPI::Allreduce(tmp.data(), globalChecks.data(), tmp.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS

        const auto invDenom = 1.0 / max(global_nPointLES[iDim], 1ul);
        mean_check_old = globalChecks[0] * invDenom;
        var_check_old = globalChecks[1] * invDenom - pow(mean_check_old, 2);

        if (!config->GetSBSParam().besselScaleFactor) {
          /*--- Standardize with the global statistics of the smoothed field (zero mean, unit variance
                over the active points), with the same clipping as the Bessel scaling. Inactive points
                stay at zero. ---*/
          const su2double invStd = 1.0/sqrt(max(var_check_old, 1e-10));
          SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
          for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
            if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
            su2double source = (nodes->GetLangevinSourceTerms(iPoint, iDim) - mean_check_old) * invStd;
            if (source < -sourceLim || source > sourceLim) source = 0.0;
            mean_check_new += source;
            var_check_new += pow(source, 2);
            nodes->SetLangevinSourceTerms(iPoint, iDim, source);
          }
          END_SU2_OMP_FOR

          /*--- Statistics after scaling, for the diagnostics. ---*/
          SU2_OMP_SAFE_GLOBAL_ACCESS(globalChecks[4] = 0.0; globalChecks[5] = 0.0;)
          atomicAdd(mean_check_new, globalChecks[4]);
          atomicAdd(var_check_new, globalChecks[5]);
          BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
            su2double tmp[2] = {globalChecks[4], globalChecks[5]};
            SU2_MPI::Allreduce(tmp, &globalChecks[4], 2, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
          }
          END_SU2_OMP_SAFE_GLOBAL_ACCESS
        }

        SU2_OMP_MASTER
        if (rank == MASTER_NODE && config->GetSBSParam().stochSourceDiagnostics) {
          mean_check_notSmoothed = globalChecks[2] * invDenom;
          var_check_notSmoothed = globalChecks[3] * invDenom - pow(mean_check_notSmoothed, 2);
          mean_check_new = globalChecks[4] * invDenom;
          var_check_new = globalChecks[5] * invDenom - pow(mean_check_new, 2);

          cout << "Mean of stochastic source term in Langevin equations:";
          cout << "\n   Uncorrelated            --> " << mean_check_notSmoothed;
          cout << "\n   Smoothed before scaling --> " << mean_check_old;
          cout << "\n   Smoothed after scaling  --> " << mean_check_new;
          cout << "\nVariance of stochastic source term in Langevin equations:";
          cout << "\n   Uncorrelated            --> " << var_check_notSmoothed;
          cout << "\n   Smoothed before scaling --> " << var_check_old;
          cout << "\n   Smoothed after scaling  --> " << var_check_new << '\n' << endl;
        }
        END_SU2_OMP_MASTER
    }
  }

}

void CTurbSASolver::SetInletAtVertex(const su2double *val_inlet,
                                    unsigned short iMarker,
                                    unsigned long iVertex) {
  SU2_ZONE_SCOPED

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];

}

su2double CTurbSASolver::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                          const CGeometry* geometry, su2double* val_inlet) const {
  SU2_ZONE_SCOPED
  const auto position = nDim + 2 + nDim;
  val_inlet[position] = Inlet_TurbVars[iMarker][iVertex][0];

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);
}

void CTurbSASolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  SU2_ZONE_SCOPED
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_TurbVars[iMarker][iVertex][0] = GetNuTilde_Inf();
    }
  }
}

void CTurbSASolver::ComputeUnderRelaxationFactor(CSolver** solver_container, const CConfig *config) {
  SU2_ZONE_SCOPED

  /* Apply the turbulent under-relaxation to the SA variants. The
   SA_NEG model is more robust due to allowing for negative nu_tilde,
   so the under-relaxation is not applied to that variant. */

  if (config->GetSAParsedOptions().version == SA_OPTIONS::NEG ||
      config->GetSBSParam().StochasticBackscatter) return;

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  const su2double allowableRatio =  config->GetMaxUpdateFractionSA();

  ComputeUnderRelaxationFactorHelper(solver_container, allowableRatio);

}
