/*!
 * \file CTurbSSTSolver.cpp
 * \brief Main subroutines of CTurbSSTSolver class
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

#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/variables/CTurbSSTVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/random_toolbox.hpp"


CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  SU2_ZONE_SCOPED
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();
  sstParsedOptions = config->GetSSTParsedOptions();

  /*--- Dimension of the problem --> dependent on the turbulence model. ---*/

  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Add Langevin equations if the Stochastic Backscatter Model is used ---*/

  if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
    nVar += 3;
    nVarGrad = nPrimVar = nVar;
  }

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0 || config->GetMGCycle() == MG_CYCLE::FULL) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy) {
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
      EdgeFluxesDiff.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
    }

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Initialize value for model constants ---*/
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  if (sstParsedOptions.version == SST_OPTIONS::V1994){
    constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
    constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
    constants[10] = 20.0; // production limiter constant
  } else {
    /* SST-V2003 */
    constants[8] = 5.0 / 9.0;  //gamma_1
    constants[9] = 0.44;  //gamma_2
    constants[10] = 10.0; // production limiter constant
  }

  /*--- Far-field flow state quantities and initialization. ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;

  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim, VelInf);

  su2double kine_Inf  = 3.0/2.0*(VelMag2*Intensity*Intensity);
  su2double omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  Solution_Inf[0] = kine_Inf;
  Solution_Inf[1] = omega_Inf;
  if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
    for (unsigned short iVar = 2; iVar < nVar; iVar++)
      Solution_Inf[iVar] = 0.0;
  }

  /*--- Constants to use for lower limit of turbulence variable. ---*/
  su2double Ck = config->GetKFactor_LowerLimit();
  su2double Cw = config->GetOmegaFactor_LowerLimit();

  /*--- Initialize lower and upper limits. ---*/
  if (sstParsedOptions.dll) {
    lowerlimit[0] = Ck * kine_Inf;
    lowerlimit[1] = Cw * omega_Inf;
  } else {
    lowerlimit[0] = 1.0e-10;
    lowerlimit[1] = 1.0e-4;
  }

  upperlimit[0] = 1.0e10;
  upperlimit[1] = 1.0e15;

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nPoint, nDim, nVar, constants, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION_EDDY);

  /*--- Initialize quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      Inlet_TurbVars[iMarker](iVertex,0) = kine_Inf;
      Inlet_TurbVars[iMarker](iVertex,1) = omega_Inf;
      if (config->GetSBSParam().StochasticBackscatter && config->GetSBSParam().stochSourceType == LANGEVIN)
        for (unsigned short iVar = 2; iVar < nVar; iVar++)
          Inlet_TurbVars[iMarker](iVertex, iVar) = 0.0;
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

  /*--- Add the solver name. ---*/
  SolverName = "SST";

}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {

    /*--- Compute the DES length scale ---*/

    SetDES_LengthScale(solver_container, geometry, config);

    bool backscatter = config->GetSBSParam().StochasticBackscatter;
    bool backscatterInBox = config->GetSBSParam().StochBackscatterInBox;
    if (backscatter && backscatterInBox) SetBackscatterInBox(config, geometry);

    /*--- Compute source terms for Langevin equations ---*/

    unsigned long innerIter = config->GetInnerIter();
    if (backscatter && innerIter==0) {
      InitiateComms(geometry, config, MPI_QUANTITIES::DES_LENGTHSCALE);
      CompleteComms(geometry, config, MPI_QUANTITIES::DES_LENGTHSCALE);
      InitiateComms(geometry, config, MPI_QUANTITIES::LES_SENSOR);
      CompleteComms(geometry, config, MPI_QUANTITIES::LES_SENSOR);
      InitiateComms(geometry, config, MPI_QUANTITIES::MEAN_TKE);
      CompleteComms(geometry, config, MPI_QUANTITIES::MEAN_TKE);

      SetLangevinSourceTerms(config, geometry);
      const unsigned short maxIter = config->GetSBSParam().SBS_maxIterSmooth;
      const su2double ctau = config->GetSBSParam().SBS_Ctau;
      if (maxIter > 0) SmoothLangevinSourceTerms(config, geometry);
      if (config->GetSBSParam().stochSourceType == ORNSTEIN_UHLENBECK) ComputeOU_Process(solver_container, config, geometry);
    }

  }

}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container,
                                    CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const su2double a1 = constants[7];

  /*--- Compute turbulence gradients. ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config, -1);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config, -1);
  }

  AD::StartNoSharedReading();

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  const su2double cYosh = 0.08;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Compute blending functions and cross diffusion ---*/

    const su2double rho = flowNodes->GetDensity(iPoint);
    const su2double mu = flowNodes->GetLaminarViscosity(iPoint);

    const su2double dist = geometry->nodes->GetWall_Distance(iPoint);

    const su2double VorticityMag = max(GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint)), 1e-12);
    const su2double StrainMag = max(flowNodes->GetStrainMag(iPoint), 1e-12);
    nodes->SetBlendingFunc(iPoint, mu, dist, rho, config->GetKind_Trans_Model());

    const su2double F2 = nodes->GetF2blending(iPoint);

    /*--- Compute the eddy viscosity ---*/

    const su2double kine = nodes->GetSolution(iPoint,0);
    const su2double omega = nodes->GetSolution(iPoint,1);

    const auto& eddy_visc_var = sstParsedOptions.version == SST_OPTIONS::V1994 ? VorticityMag : StrainMag;
    const su2double muT = max(0.0, rho * a1 * kine / max(a1 * omega, eddy_visc_var * F2));

    nodes->SetmuT(iPoint, muT);

  }
  END_SU2_OMP_FOR


  /*--- Compute turbulence index ---*/
  if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
    for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (!config->GetViscous_Wall(iMarker)) continue;

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

        if (!geometry->nodes->GetDomain(iPoint)) continue;

        const auto jPoint = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        su2double shearStress = 0.0;
        for(auto iDim = 0u; iDim < nDim; iDim++) {
          shearStress += pow(solver_container[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, iDim), 2.0);
        }
        shearStress = sqrt(shearStress);

        const su2double FrictionVelocity = max(sqrt(shearStress/flowNodes->GetDensity(iPoint)), EPS);
        const su2double wall_dist = geometry->vertex[iMarker][iVertex]->GetNearestNeighborDistance();

        const su2double Derivative = flowNodes->GetLaminarViscosity(jPoint) * pow(nodes->GetSolution(jPoint, 0), 0.673) / wall_dist;
        const su2double turbulence_index = 6.1 * Derivative / pow(FrictionVelocity, 2.346);

        nodes->SetTurbIndex(iPoint, turbulence_index);
      }
      END_SU2_OMP_FOR
    }
  }

  AD::EndNoSharedReading();
}

void CTurbSSTSolver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, const CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/
  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {
    /*--- Menter's first blending function (only SST)---*/
    numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(jPoint));
  };

  /*--- Now instantiate the generic non-conservative implementation with the functor above. ---*/
  Viscous_Residual_NonCons(iEdge, geometry, solver_container, numerics, config, SolverSpecificNumerics);

}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  bool axisymmetric = config->GetAxisymmetric();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Menter's first blending function ---*/

    numerics->SetF1blending(nodes->GetF1blending(iPoint),0.0);

    /*--- Menter's second blending function ---*/

    numerics->SetF2blending(nodes->GetF2blending(iPoint));

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Cross diffusion ---*/

    numerics->SetCrossDiff(nodes->GetCrossDiff(iPoint));

    /*--- Effective Intermittency ---*/
    if (config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      numerics->SetIntermittencyEff(solver_container[TRANS_SOL]->GetNodes()->GetIntermittencyEff(iPoint));
    }

    if (axisymmetric){
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(iPoint));
    }

    if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
      su2double k = nodes->GetSolution(iPoint, 0);
      su2double omega = nodes->GetSolution(iPoint, 1);
      su2double beta_star = constants[6];
      su2double RANS_lengthscale = sqrt(k)/(beta_star*max(omega, 1e-10));
      su2double DES_lengthscale = nodes->GetDES_LengthScale(iPoint);
      su2double FDDES = RANS_lengthscale/max(DES_lengthscale, 1e-10);
      numerics->SetFDDES(FDDES, 0.0);

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
        numerics->SetMaxDelta(geometry->nodes->GetMaxLength(iPoint), 0.0);
        if (config->GetSBSParam().useMeanTurbKE) numerics->SetAvgTurbKineticEnergy(nodes->GetMeanTurbKinEnergy(iPoint), 0.0);
      }
    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Store the intermittency ---*/

    if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
      nodes->SetIntermittency(iPoint, numerics->GetIntermittencyEff());
    }

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  /*--- Custom user defined source term (from the python wrapper) ---*/
  if (config->GetPyCustomSource()) {
    CustomSourceResidual(geometry, solver_container, numerics_container, config, iMesh);
  }

}

void CTurbSSTSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  WALL_TYPE WallType; su2double Roughness_Height;
  tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
  const bool rough_wall = WallType == WALL_TYPE::ROUGH;

  /*--- Evaluate nu tilde at the closest point to the surface using the wall functions. ---*/

  if (config->GetWall_Functions()) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(SetTurbVars_WF(geometry, solver_container, config, val_marker);)
    return;
  }

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- distance to closest neighbor ---*/
    su2double wall_dist = geometry->vertex[val_marker][iVertex]->GetNearestNeighborDistance();

    su2double solution[MAXNVAR];

    if (rough_wall) {
      /*--- Set wall values ---*/
      su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
      su2double WallShearStress = solver_container[FLOW_SOL]->GetWallShearStress(val_marker, iVertex);

      /*--- Compute non-dimensional velocity ---*/
      su2double FrictionVel = sqrt(fabs(WallShearStress)/density);

      /*--- Compute roughness in wall units. ---*/
      su2double kPlus = FrictionVel*Roughness_Height*density/laminar_viscosity;

      /*--- Modify the omega and k to account for a rough wall. ---*/

      switch (config->GetKindRoughSSTModel()) {
        /*--- Reference 1 original Wilcox (1998). ---*/
        case ROUGHSST_MODEL::WILCOX1998: {
          su2double S_R = 0.0;
          if (kPlus <= 25)
            S_R = pow(50/(kPlus+EPS), 2);
          else
            S_R = 100/(kPlus+EPS);

          solution[0] = 0.0;
          solution[1] = FrictionVel*FrictionVel*S_R/(laminar_viscosity/density);
        } break;
        /*--- Reference 2 from D.C. Wilcox Turbulence Modeling for CFD (2006) ---*/
        case ROUGHSST_MODEL::WILCOX2006: {
          su2double S_R = 0.0;
          if (kPlus <= 5)
            S_R = pow(200/(kPlus+EPS),2);
          else
            S_R = 100/(kPlus+EPS) + (pow(200/(kPlus+EPS),2) - 100/(kPlus+EPS))*exp(5-kPlus);

          solution[0] = 0.0;
          solution[1] = FrictionVel*FrictionVel*S_R/(laminar_viscosity/density);
        } break;
        /*--- Knopp eddy viscosity limiter ---*/
        case ROUGHSST_MODEL::LIMITER_KNOPP: {
          su2double d0 = 0.03*Roughness_Height*min(1.0, pow((kPlus + EPS )/30.0, 2.0/3.0))*min(1.0, pow((kPlus + EPS)/45.0, 0.25))*min(1.0, pow((kPlus + EPS) /60, 0.25));
          solution[0] = (FrictionVel*FrictionVel / sqrt(constants[6]))*min(1.0, kPlus / 90.0);

          const su2double kappa = config->GetwallModel_Kappa();
          su2double beta_1 = constants[4];
          solution[1] = min( FrictionVel/(sqrt(constants[6])*d0*kappa), 60.0*laminar_viscosity/(density*beta_1*pow(wall_dist,2)));
        } break;
        /*--- Aupoix eddy viscosity limiter ---*/
        case (ROUGHSST_MODEL::LIMITER_AUPOIX): {
          su2double k0Plus = ( 1.0 /sqrt( constants[6])) * tanh((log10((kPlus +EPS ) / 30.0) + 1.0 - 1.0*tanh( (kPlus + EPS) / 125.0))*tanh((kPlus + EPS) / 125.0));
          su2double kwallPlus = max(0.0, k0Plus);
          su2double kwall = kwallPlus*FrictionVel*FrictionVel;

          su2double omegawallPlus = (300.0 / pow(kPlus + EPS, 2.0)) * pow(tanh(15.0 / (4.0*kPlus)), -1.0) + (191.0 / (kPlus + EPS))*(1.0 - exp(-kPlus / 250.0));

          solution[0] = kwall;
          solution[1] = omegawallPlus*FrictionVel*FrictionVel*density/laminar_viscosity;
        } break;
      }
    } else { // smooth wall
      /*--- Set wall values ---*/
      su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      su2double laminar_viscosity = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

      su2double beta_1 = constants[4];
      solution[0] = 0.0;
      solution[1] = 60.0*laminar_viscosity/(density*beta_1*pow(wall_dist,2));
    }

    /*--- Set the solution values and zero the residual ---*/
    nodes->SetSolution_Old(iPoint, solution);
    nodes->SetSolution(iPoint, solution);
    LinSysRes.SetBlock_Zero(iPoint);

    if (implicit) {
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint, 0);
      Jacobian.DeleteValsRowi(iPoint, 1);
    }
  }
  END_SU2_OMP_FOR
}


void CTurbSSTSolver::SetTurbVars_WF(CGeometry *geometry, CSolver **solver_container,
                                    const CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- von Karman constant from boundary layer theory ---*/
  const su2double kappa = config->GetwallModel_Kappa();
  const su2double minYPlus = config->GetwallModel_MinYPlus();
  /*--- relaxation factor for k-omega values ---*/
  const su2double relax = config->GetwallModel_RelFac();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    const auto iPoint_Neighbor = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    if (!geometry->nodes->GetDomain(iPoint_Neighbor)) continue;

    su2double Y_Plus = solver_container[FLOW_SOL]->GetYPlus(val_marker, iVertex);
    su2double Lam_Visc_Wall = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

    /*--- Do not use wall model at the ipoint when y+ < "limit", use zero flux (Neumann) conditions. ---*/

    if (Y_Plus < minYPlus) {
      /* --- Use zero flux (Neumann) conditions, i.e. nothing has to be done. --- */
      continue;
    }

    su2double Eddy_Visc = solver_container[FLOW_SOL]->GetEddyViscWall(val_marker, iVertex);
    su2double k = nodes->GetSolution(iPoint_Neighbor,0);
    su2double omega = nodes->GetSolution(iPoint_Neighbor,1);
    su2double Density_Wall = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    su2double U_Tau = solver_container[FLOW_SOL]->GetUTau(val_marker, iVertex);
    su2double y = Y_Plus*Lam_Visc_Wall/(Density_Wall*U_Tau);

    su2double omega1 = 6.0*Lam_Visc_Wall/(0.075*Density_Wall*y*y);  // eq. 19
    su2double omega0 = U_Tau/(sqrt(0.09)*kappa*y);                  // eq. 20
    su2double omega_new = sqrt(omega0*omega0 + omega1*omega1);      // eq. 21 Nichols & Nelson
    su2double k_new = omega_new * Eddy_Visc/Density_Wall;           // eq. 22 Nichols & Nelson
                                           // (is this the correct density? paper says rho and not rho_w)

    /*--- put some relaxation factor on the k-omega values ---*/
    k += relax*(k_new - k);
    omega += relax*(omega_new - omega);

    su2double solution[MAXNVAR] = {k, omega};

    nodes->SetSolution_Old(iPoint_Neighbor,solution);
    nodes->SetSolution(iPoint,solution);

    LinSysRes.SetBlock_Zero(iPoint_Neighbor);

    if (implicit) {
      /*--- includes 1 in the diagonal ---*/
      Jacobian.DeleteValsRowi(iPoint_Neighbor, 0);
      Jacobian.DeleteValsRowi(iPoint_Neighbor, 1);
    }
  }
}

void CTurbSSTSolver::SetDES_LengthScale(CSolver **solver, CGeometry *geometry, CConfig *config){
  SU2_ZONE_SCOPED

  const su2double cDES_1 = config->GetConst_DES_1();
  const su2double cDES_2 = config->GetConst_DES_2();
  const su2double beta_star = constants[6];
  const su2double k2 = pow(0.41, 2);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver[FLOW_SOL]->GetNodes());

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++){

    const auto wallDistance  = geometry->nodes->GetWall_Distance(iPoint);
    const auto velocityGrad  = flowNodes->GetVelocityGradient(iPoint);
    const auto density       = flowNodes->GetDensity(iPoint);
    const auto laminarViscosity = flowNodes->GetLaminarViscosity(iPoint);
    const auto eddyViscosity    = nodes->GetmuT(iPoint);
    const su2double kinematicViscosity     = laminarViscosity/density;
    const su2double kinematicViscosityTurb = eddyViscosity/density;
    const su2double k = nodes->GetSolution(iPoint, 0);
    const su2double omega = max(nodes->GetSolution(iPoint, 1), 1e-10);
    const su2double F1 = nodes->GetF1blending(iPoint);
    su2double constDES = cDES_1*F1 + cDES_2*(1.0-F1);

    su2double uijuij = 0.0;
    for(auto iDim = 0u; iDim < nDim; iDim++){
      for(auto jDim = 0u; jDim < nDim; jDim++){
        uijuij += pow(velocityGrad[iDim][jDim], 2);
      }
    }
    uijuij = sqrt(fabs(uijuij));
    uijuij = max(uijuij,1e-10);

    const su2double LES_FilterWidth = config->GetLES_FilterWidth();
    su2double maxDelta = (LES_FilterWidth > 0.0) ? LES_FilterWidth : geometry->nodes->GetMaxLength(iPoint);

    const su2double r_d = (kinematicViscosityTurb+kinematicViscosity)/(uijuij*k2*pow(wallDistance, 2.0));
    const su2double f_d = 1.0-tanh(pow(20.0*r_d,3.0));

    const su2double distDES = constDES * maxDelta;
    const su2double distRANS = sqrt(k)/(omega*beta_star);
    su2double lengthScale = distRANS-f_d*max(0.0,(distRANS-distDES));
    su2double lesSensor = (distRANS<=distDES) ? 0.0 : f_d;

    if (config->GetEnforceLES()) {
      lengthScale = distDES;
      lesSensor = 1.0;
    }

    nodes->SetDES_LengthScale(iPoint, lengthScale);
    nodes->SetLES_Mode(iPoint, lesSensor);

  }
  END_SU2_OMP_FOR
}

void CTurbSSTSolver::SetBackscatterInBox(CConfig *config, CGeometry *geometry) {
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

void CTurbSSTSolver::SetLangevinSourceTerms(CConfig *config, CGeometry* geometry) {
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
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          nodes->SetLangevinSourceTermsOld(iPoint, iDim, dummySource);
          nodes->SetLangevinSourceTerms(iPoint, iDim, 0.0);
        }
      }
    }
    END_SU2_OMP_FOR
  }
}

void CTurbSSTSolver::ComputeOU_Process(CSolver **solver, CConfig *config, CGeometry *geometry) {
  SU2_ZONE_SCOPED

  su2double timeStep = config->GetDelta_UnstTimeND();
  const su2double beta = constants[6];

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double maxDelta = geometry->nodes->GetMaxLength(iPoint);
    su2double tke = (config->GetSBSParam().useMeanTurbKE) ? nodes->GetMeanTurbKinEnergy(iPoint) : nodes->GetSolution(iPoint, 0);
    su2double timeLES = fabs(config->GetSBSParam().SBS_Ctau) * maxDelta / sqrt(max(tke, 1e-10));
    su2double timeRANS = 1.0 / (beta*max(nodes->GetSolution(iPoint, 1), 1e-10));
    su2double lesMode = nodes->GetLES_Mode(iPoint);
    su2double timeBlended = min(timeRANS, 10.0*timeStep) * (1.0-lesMode) + timeLES * lesMode;
    su2double timeRatio = timeStep/timeBlended;
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

  InitiateComms(geometry, config, MPI_QUANTITIES::OU_PROCESS);
  CompleteComms(geometry, config, MPI_QUANTITIES::OU_PROCESS);
}

void CTurbSSTSolver::SmoothLangevinSourceTerms(CConfig* config, CGeometry* geometry) {
  SU2_ZONE_SCOPED

  static su2double globalResNorm;
  static unsigned long global_nPointLES;
  static std::array<su2double, 6> globalChecks;
  /*--- Scalars driving BiCGSTAB control flow (loop conditions/break). These must be "static" (as
        in the original Jacobi implementation) rather than plain locals: the whole point-loop and
        Allreduce machinery below runs master-thread-only per MPI rank (BEGIN/END_SU2_OMP_SAFE_
        GLOBAL_ACCESS), but the while-loop conditions and breakdown checks are plain code executed
        redundantly by every OpenMP thread sharing this parallel region, so they must all observe
        the exact same, already-synchronized value that only the master thread computed. ---*/
  static su2double globalRho, globalRhoPrev, globalAlpha, globalOmega, globalR0V, globalTS, globalTT;
  static bool breakdown;

  const su2double LES_FilterWidth = config->GetLES_FilterWidth();
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
        diffusive flux across each face (deferred correction, folded into the RHS, see below). */

  if (timeIter == restartIter) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      su2double maxDelta = (LES_FilterWidth > 0.0) ? LES_FilterWidth : geometry->nodes->GetMaxLength(iPoint);
      su2double b2 = cDelta * maxDelta * maxDelta;
      su2double volume_iPoint = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      auto coord_i = geometry->nodes->GetCoord(iPoint);
      su2double diag = 1.0;
      for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
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
        diag += a_ij;

        /*--- betaVec = (b2/V_i) * sign * (normal - (dot_nd/dist_ij_2)*edge_vector), the coefficient
              such that mean_grad_face . betaVec gives the non-orthogonal (tangential) part of the
              gradient-based diffusive flux, consistent with CAvgGrad_Base::CorrectGradient. ---*/

        for (unsigned short index = 0; index < nDim; index++) {
          su2double betaVec = (b2/volume_iPoint) * sign * (normal[index] - (dot_nd/dist_ij_2) * dx_ij_vec[index]);
          nodes->SetSmoothingBetaVec(iPoint, iNode, index, betaVec);
        }
      }
      nodes->SetSmoothingDiag(iPoint, diag);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Matrix-free operator A(x)_i = diag_i*x_i - sum_j(a_ij*x_j), the only part of the discrete
        system that is a genuine matrix (the tangential/non-orthogonal term below is a fixed RHS
        contribution between gradient refreshes, not part of A). GetInput reads whichever field is
        currently being multiplied by A (the solution itself to form the initial residual, or the
        preconditioned BiCGSTAB directions phat/shat); the result is written into a local, rank-
        private (never used as a neighbor, so never halo-exchanged) work vector. Runs master-thread-
        only per rank, exactly like the original Jacobi sweep. ---*/

  auto ComputeMatVec = [&](unsigned short iDim, auto&& GetInput, std::vector<su2double>& out) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
      su2double diag = nodes->GetSmoothingDiag(iPoint);
      su2double sum = 0.0;
      for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
        auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
        su2double a_ij = nodes->GetSmoothingMatrixCoeff(iPoint, iNode);
        sum += a_ij * GetInput(jPoint);
      }
      out[iPoint] = diag*GetInput(iPoint) - sum;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  };

  /*--- Dot product between two rank-private work vectors, restricted to the points actually being
        solved for (points outside the active LES region are frozen at 0 and never enter the linear
        system, mirroring the exclusion in the original Jacobi sweep), reduced over all MPI ranks.
        Writes into the static scalar "out" so every OpenMP thread sharing this parallel region
        observes the same, already-synchronized value once outside the master-only block. ---*/

  auto DotActive = [&](unsigned short iDim, const std::vector<su2double>& a, const std::vector<su2double>& b,
                       su2double& out) {
    su2double local = 0.0;
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
      local += a[iPoint]*b[iPoint];
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    SU2_MPI::Allreduce(&local, &out, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  };

  /*--- Solve, for each spatial dimension, diag_i*x_i - sum_j(a_ij*x_j) = source_i_old + tangential_i
        with a matrix-free, Jacobi(diagonal)-preconditioned BiCGSTAB, restarting the Krylov subspace
        (warm-started from the current solution) every time the tangential/gradient RHS term is
        refreshed, so the matrix and RHS are always fixed for the duration of a Krylov subspace. ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {

    std::vector<su2double> bRhs(nPointDomain, 0.0), rVec(nPointDomain, 0.0), rhat0Vec(nPointDomain, 0.0);
    std::vector<su2double> pVec(nPointDomain, 0.0), vVec(nPointDomain, 0.0);
    std::vector<su2double> sVec(nPointDomain, 0.0), tVec(nPointDomain, 0.0);

    unsigned short totalIter = 0;
    bool converged = false;

    SU2_OMP_MASTER
    if (rank == MASTER_NODE) {
      cout << "\nResidual of Laplacian smoothing along dimension " << iDim+1
           << "\n---------------------------------"
           << "\n   Iter       RMS Residual"
           << "\n---------------------------------" << endl;
    }
    END_SU2_OMP_MASTER

    while (!converged && totalIter < maxIter) {

      /*--- Start (or restart) of a block: refresh the halo values of the solution, recompute the
            Green-Gauss gradient used for the non-orthogonal RHS correction, and rebuild the initial
            residual r0 = b - A*x (warm-started from the current solution estimate). ---*/

      InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);
      CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG);

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        su2double source_i = nodes->GetLangevinSourceTerms(iPoint, iDim);
        su2double volume_iPoint = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
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
          nodes->SetLangevinSourceGrad(iPoint, index, grad_i[index] / max(volume_iPoint, 1e-10));
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      InitiateComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG_GRAD);
      CompleteComms(geometry, config, MPI_QUANTITIES::STOCH_SOURCE_LANG_GRAD);

      unsigned long local_nPointLES = 0;

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        su2double source_i_old = nodes->GetLangevinSourceTermsOld(iPoint, iDim);
        if (source_i_old > 3.0*sourceLim) continue;
        local_nPointLES += 1;

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
        bRhs[iPoint] = source_i_old + tangential_i;
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      SU2_MPI::Allreduce(&local_nPointLES, &global_nPointLES, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      ComputeMatVec(iDim, [&](unsigned long j){ return nodes->GetLangevinSourceTerms(j, iDim); }, rVec);

      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
        rVec[iPoint] = bRhs[iPoint] - rVec[iPoint];
        rhat0Vec[iPoint] = rVec[iPoint];
        /*--- p_0 = v_0 = 0 (standard BiCGSTAB init): with rho_prev=alpha=omega=1 reset below, the
              first inner iteration's p update (p = r + beta*(p-omega*v)) then correctly reduces to
              p_1 = r_0 regardless of beta, since the (p_0 - omega*v_0) term vanishes. ---*/
        pVec[iPoint] = 0.0;
        vVec[iPoint] = 0.0;
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS

      /*--- Account for the cost of this refresh/restart against the iteration budget: this
            guarantees the outer (block/restart) loop always makes progress towards maxIter, even
            in the degenerate case where every restarted Krylov subspace breaks down immediately
            (e.g. an already-converged residual, for which rho would be ~0 by construction). ---*/
      totalIter++;

      su2double r0NormSq = 0.0;
      DotActive(iDim, rVec, rVec, r0NormSq);
      SU2_OMP_SAFE_GLOBAL_ACCESS(
        globalResNorm = (global_nPointLES==0) ? su2double(0.0) : sqrt(r0NormSq / global_nPointLES);
      )
      if (log10(globalResNorm) < tol) converged = true;

      SU2_OMP_SAFE_GLOBAL_ACCESS(globalRho = 1.0; globalAlpha = 1.0; globalOmega = 1.0; breakdown = false;)

      unsigned short blockIter = 0;
      while (!breakdown && !converged && blockIter < gradRefreshInterval && totalIter < maxIter) {

        SU2_OMP_SAFE_GLOBAL_ACCESS(globalRhoPrev = globalRho;)
        DotActive(iDim, rhat0Vec, rVec, globalRho);

        if (fabs(globalRho) < eps_breakdown) { breakdown = true; break; }

        su2double beta = (globalRho/globalRhoPrev) * (globalAlpha/globalOmega);

        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
          pVec[iPoint] = rVec[iPoint] + beta*(pVec[iPoint] - globalOmega*vVec[iPoint]);
          nodes->SetSmoothPhat(iPoint, iDim, pVec[iPoint] / nodes->GetSmoothingDiag(iPoint));
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS

        InitiateComms(geometry, config, MPI_QUANTITIES::SMOOTH_PHAT);
        CompleteComms(geometry, config, MPI_QUANTITIES::SMOOTH_PHAT);

        ComputeMatVec(iDim, [&](unsigned long j){ return nodes->GetSmoothPhat(j, iDim); }, vVec);

        DotActive(iDim, rhat0Vec, vVec, globalR0V);
        if (fabs(globalR0V) < eps_breakdown) { breakdown = true; break; }
        SU2_OMP_SAFE_GLOBAL_ACCESS(globalAlpha = globalRho / globalR0V;)

        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
          sVec[iPoint] = rVec[iPoint] - globalAlpha*vVec[iPoint];
          nodes->SetSmoothShat(iPoint, iDim, sVec[iPoint] / nodes->GetSmoothingDiag(iPoint));
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS

        InitiateComms(geometry, config, MPI_QUANTITIES::SMOOTH_SHAT);
        CompleteComms(geometry, config, MPI_QUANTITIES::SMOOTH_SHAT);

        ComputeMatVec(iDim, [&](unsigned long j){ return nodes->GetSmoothShat(j, iDim); }, tVec);

        DotActive(iDim, tVec, sVec, globalTS);
        DotActive(iDim, tVec, tVec, globalTT);
        if (fabs(globalTT) < eps_breakdown) { breakdown = true; break; }
        SU2_OMP_SAFE_GLOBAL_ACCESS(globalOmega = globalTS / globalTT;)

        su2double localResNorm = 0.0;
        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
          if (nodes->GetLangevinSourceTermsOld(iPoint, iDim) > 3.0*sourceLim) continue;
          su2double source_i = nodes->GetLangevinSourceTerms(iPoint, iDim);
          su2double phat = nodes->GetSmoothPhat(iPoint, iDim);
          su2double shat = nodes->GetSmoothShat(iPoint, iDim);
          source_i += globalAlpha*phat + globalOmega*shat;
          nodes->SetLangevinSourceTerms(iPoint, iDim, source_i);
          rVec[iPoint] = sVec[iPoint] - globalOmega*tVec[iPoint];
          localResNorm += pow(rVec[iPoint], 2);
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS

        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
          SU2_MPI::Allreduce(&localResNorm, &globalResNorm, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
          globalResNorm = (global_nPointLES==0) ? su2double(0.0) : sqrt(globalResNorm / global_nPointLES);
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS

        blockIter++; totalIter++;

        SU2_OMP_MASTER
        if (rank == MASTER_NODE) {
          cout << "  "
               << std::setw(5) << totalIter-1
               << "       "
               << std::setw(12) << std::fixed << std::setprecision(6) << log10(globalResNorm)
               << endl;
        }
        END_SU2_OMP_MASTER

        if (log10(globalResNorm) < tol || totalIter == maxIter) converged = true;
      }
    }

    SU2_OMP_MASTER
    if (rank == MASTER_NODE) cout << "---------------------------------" << endl;
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
            su2double integral = 0.0;
            if (timeIter==restartIter) {
              su2double maxDelta = (LES_FilterWidth > 0.0) ? LES_FilterWidth : geometry->nodes->GetMaxLength(iPoint);
              su2double b2 = cDelta * maxDelta * maxDelta;
              su2double M[3][3] = {{0.0}};
              for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
                auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
                auto iEdge = geometry->nodes->GetEdge(iPoint, iNode);
                auto* normal = geometry->edges->GetNormal(iEdge);
                for (unsigned short ind1 = 0; ind1 < nDim; ind1++) {
                  for (unsigned short ind2 = 0; ind2 < nDim; ind2++) {
                    M[ind1][ind2] += normal[ind1]*normal[ind2];
                  }
                }
              }
              su2double a = M[0][0], b = M[1][1], c = M[2][2];
              su2double d = M[0][1], e = M[1][2], f = M[0][2];
              su2double lambda[3] = {0.0};
              su2double p1 = d*d + e*e + f*f;
              if (p1 < 1e-20) {
                lambda[0] = a;
                lambda[1] = b;
                lambda[2] = (nDim==3 ? c : b);
              } else {
                su2double trace = (a + b + c) / 3.0;
                su2double p2 = (a-trace)*(a-trace) +
                               (b-trace)*(b-trace) +
                               (c-trace)*(c-trace) +
                               2.0 * p1;
                su2double p = sqrt(p2 / 6.0);
                su2double B[3][3];
                for (unsigned short ind1 = 0; ind1 < nDim; ind1++)
                  for (unsigned short ind2 = 0; ind2 < nDim; ind2++)
                    B[ind1][ind2] = M[ind1][ind2];
                B[0][0] -= trace;
                B[1][1] -= trace;
                B[2][2] -= trace;
                for (unsigned short ind1 = 0; ind1 < nDim; ind1++)
                  for (unsigned short ind2 = 0; ind2 < nDim; ind2++)
                    B[ind1][ind2] /= p;
                su2double detB =
                    B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[2][1]) -
                    B[0][1]*(B[1][0]*B[2][2] - B[1][2]*B[2][0]) +
                    B[0][2]*(B[1][0]*B[2][1] - B[1][1]*B[2][0]);
                su2double r = detB * 0.5;
                r = max(min(r, 1.0), -1.0);
                su2double phi = acos(r) / 3.0;
                lambda[0] = max(trace + 2.0*p*cos(phi), 1e-10);
                lambda[1] = max(trace + 2.0*p*cos(phi + 2.0*M_PI/3.0), 1e-10);
                lambda[2] = max(trace + 2.0*p*cos(phi + 4.0*M_PI/3.0), 1e-10);
              }
              su2double V = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
              su2double dI = V / sqrt(lambda[0]);
              su2double dJ = V / sqrt(lambda[1]);
              su2double dK = V / sqrt(lambda[2]);
              su2double dI2 = dI * dI;
              su2double dJ2 = dJ * dJ;
              su2double dK2 = dK * dK;
              su2double bI = b2 / dI2;
              su2double bJ = b2 / dJ2;
              su2double bK = b2 / dK2;
              integral = RandomToolbox::GetBesselIntegral(bI, bJ, bK);
              nodes->SetBesselIntegral(iPoint, integral);
            } else {
              integral = nodes->GetBesselIntegral(iPoint);
            }
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

        const auto invDenom = 1.0 / max(global_nPointLES, 1ul);
        mean_check_old = globalChecks[0] * invDenom;
        var_check_old = globalChecks[1] * invDenom - pow(mean_check_old, 2);

        if (!config->GetSBSParam().besselScaleFactor) {
          SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
          for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
            su2double source = nodes->GetLangevinSourceTerms(iPoint, iDim);
            source *= 1.0/sqrt(max(var_check_old, 1e-10));
            nodes->SetLangevinSourceTerms(iPoint, iDim, source);
          }
          END_SU2_OMP_FOR
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
          cout << "\n   Smoothed after scaling  --> " << ((config->GetSBSParam().besselScaleFactor) ? var_check_new : 1.0) << '\n' << endl;
        }
        END_SU2_OMP_MASTER
    }
  }

}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                              CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      su2double Inlet_Vars[MAXNVAR];
      if (config->GetInlet_Profile_From_File()) {
        /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
        Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0] / pow(config->GetVelocity_Ref(), 2);
        Inlet_Vars[1] = Inlet_TurbVars[val_marker][iVertex][1] * config->GetViscosity_Ref() /
                        (config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2));
      } else {
        /*--- Obtain fluid model for computing the  kine and omega to impose at the inlet boundary. ---*/
        CFluidModel* FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

        /*--- Obtain flow velocity vector at inlet boundary node ---*/

        const su2double* Velocity_Inlet = &V_inlet[prim_idx.Velocity()];
        su2double Density_Inlet;
        if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
          Density_Inlet = V_inlet[prim_idx.Density()];
          FluidModel->SetTDState_Prho(V_inlet[prim_idx.Pressure()], Density_Inlet);
        } else {
          const su2double* Scalar_Inlet = nullptr;
          if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
            Scalar_Inlet = config->GetInlet_SpeciesVal(config->GetMarker_All_TagBound(val_marker));
          }
          FluidModel->SetTDState_T(V_inlet[prim_idx.Temperature()], Scalar_Inlet);
          Density_Inlet = FluidModel->GetDensity();
        }
        const su2double Laminar_Viscosity_Inlet = FluidModel->GetLaminarViscosity();
        const su2double* Turb_Properties = config->GetInlet_TurbVal(config->GetMarker_All_TagBound(val_marker));
        const su2double Intensity = Turb_Properties[0];
        const su2double viscRatio = Turb_Properties[1];
        const su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim, Velocity_Inlet);

        Inlet_Vars[0] = 3.0 / 2.0 * (VelMag2 * pow(Intensity, 2));
        Inlet_Vars[1] = Density_Inlet * Inlet_Vars[0] / (Laminar_Viscosity_Inlet * viscRatio);
      }

      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/
      /*--- Load the inlet turbulence variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_inlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      //
      //      su2double Coord_Reflected[MAXNDIM];
      //      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
      //                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      //      visc_numerics->SetNormal(Normal);
      //
      //      /*--- Conservative variables w/o reconstruction ---*/
      //
      //      visc_numerics->SetPrimitive(V_domain, V_inlet);
      //
      //      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      //
      //     visc_numerics->SetScalarVar(Solution_i, Solution_j);
      //     visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      //
      //      /*--- Menter's first blending function ---*/
      //
      //      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
      //
      //      /*--- Compute residual, and Jacobians ---*/
      //
      //      auto residual = visc_numerics->ComputeResidual(config);
      //
      //      /*--- Subtract residual, and update Jacobians ---*/
      //
      //      LinSysRes.SubtractBlock(iPoint, residual);
      //      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }

  }
  END_SU2_OMP_FOR
}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the outlet ---*/

      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set Normal (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_outlet[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems ---*/
//
//      su2double Coord_Reflected[MAXNDIM];
//      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//
//      visc_numerics->SetScalarVar(Solution_i, Solution_j);
//      visc_numerics->SetScalarVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Menter's first blending function ---*/
//
//      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(), node[iPoint]->GetF1blending());
//
//      /*--- Compute residual, and Jacobians ---*/
//
//      auto residual = visc_numerics->ComputeResidual(config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//
//      LinSysRes.SubtractBlock(iPoint, residual);
//      Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR
}


void CTurbSSTSolver::BC_Inlet_MixingPlane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                          CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    su2double extAverageKine = solver_container[FLOW_SOL]->GetExtAverageKine(val_marker, iSpan);
    su2double extAverageOmega = solver_container[FLOW_SOL]->GetExtAverageOmega(val_marker, iSpan);
    su2double solution_j[] = {extAverageKine, extAverageOmega};

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/
      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSTSolver::BC_Inlet_Turbo(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  const auto nSpanWiseSections = config->GetnSpanWiseSections();

  /*--- Quantities for computing the  kine and omega to impose at the inlet boundary. ---*/
  CFluidModel *FluidModel = solver_container[FLOW_SOL]->GetFluidModel();

  su2double Intensity = config->GetTurbulenceIntensity_FreeStream();
  su2double viscRatio = config->GetTurb2LamViscRatio_FreeStream();

  for (auto iSpan = 0u; iSpan < nSpanWiseSections ; iSpan++){

    /*--- Compute the inflow kine and omega using the span wise averge quntities---*/

    su2double rho       = solver_container[FLOW_SOL]->GetAverageDensity(val_marker, iSpan);
    su2double pressure  = solver_container[FLOW_SOL]->GetAveragePressure(val_marker, iSpan);
    su2double kine      = solver_container[FLOW_SOL]->GetAverageKine(val_marker, iSpan);

    FluidModel->SetTDState_Prho(pressure, rho);
    su2double muLam = FluidModel->GetLaminarViscosity();

    su2double VelMag2 = GeometryToolbox::SquaredNorm(nDim,
                          solver_container[FLOW_SOL]->GetAverageTurboVelocity(val_marker, iSpan));

    su2double kine_b  = 3.0/2.0*(VelMag2*Intensity*Intensity);
    su2double omega_b = rho*kine/(muLam*viscRatio);

    su2double solution_j[] = {kine_b, omega_b};

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {

      /*--- find the node related to the vertex ---*/
      const auto iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();

      /*--- using the other vertex information for retrieving some information ---*/
      const auto oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();

      /*--- Index of the closest interior node ---*/
      const auto Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][oldVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, oldVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states. Use average span-wise values
             values for the turbulent state at the inflow. ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto conv_residual = conv_numerics->ComputeResidual(config);

      /*--- Jacobian contribution for implicit integration ---*/
      LinSysRes.AddBlock(iPoint, conv_residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, conv_residual.jacobian_i);

      /*--- Viscous contribution ---*/
      su2double Coord_Reflected[MAXNDIM];
      GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                               geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(nodes->GetF1blending(iPoint), nodes->GetF1blending(iPoint));

      /*--- Compute residual, and Jacobians ---*/
      auto visc_residual = visc_numerics->ComputeResidual(config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, visc_residual);
      if (implicit) Jacobian.SubtractBlock2Diag(iPoint, visc_residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

}

void CTurbSSTSolver::SetInletAtVertex(const su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {
  SU2_ZONE_SCOPED

  Inlet_TurbVars[iMarker][iVertex][0] = val_inlet[nDim+2+nDim];
  Inlet_TurbVars[iMarker][iVertex][1] = val_inlet[nDim+2+nDim+1];

}

su2double CTurbSSTSolver::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                           const CGeometry* geometry, su2double* val_inlet) const {
  SU2_ZONE_SCOPED
  const auto tke_position = nDim + 2 + nDim;
  const auto omega_position = tke_position + 1;
  val_inlet[tke_position] = Inlet_TurbVars[iMarker][iVertex][0];
  val_inlet[omega_position] = Inlet_TurbVars[iMarker][iVertex][1];

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);}

void CTurbSSTSolver::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  SU2_ZONE_SCOPED
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_TurbVars[iMarker][iVertex][0] = GetTke_Inf();
      Inlet_TurbVars[iMarker][iVertex][1] = GetOmega_Inf();
    }
  }

}

void CTurbSSTSolver::ComputeUnderRelaxationFactor(CSolver** solver_container, const CConfig *config) {
  SU2_ZONE_SCOPED

  if (config->GetSBSParam().StochasticBackscatter) return;

  const su2double allowableRatio = config->GetMaxUpdateFractionSST();

  ComputeUnderRelaxationFactorHelper(solver_container, allowableRatio);
}