/*!
 * \file CFlowOutput.hpp
 * \brief  Headers of the flow output.
 * \author F. Palacios, T. Economon, M. Colonno
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

#pragma once

#include "CFVMOutput.hpp"
#include "../variables/CVariable.hpp"

/*--- Forward declare to avoid including here. ---*/
template <class>
struct CPrimitiveIndices;

class CFlowOutput : public CFVMOutput{
protected:
  unsigned long lastInnerIter;

  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowOutput(const CConfig *config, unsigned short nDim, bool femOutput);

  /*!
   * \brief Add turboperformance outputs as history field
   * \param[in] nZone - Number of zones in problem
  */
  void AddTurboOutput(unsigned short nZone);

  /*!
   * \brief Set the values of the volume output fields for a surface point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iPoint - Index of the point.
   * \param[in] iMarker - Index of the surface marker.
   * \param[in] iVertex - Index of the vertex on the marker.
   */
  void LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver,
                       unsigned long iPoint, unsigned short iMarker, unsigned long iVertex) override;

  /*!
   * \brief Add flow surface output fields
   * \param[in] config - Definition of the particular problem.
   */
  void AddAnalyzeSurfaceOutput(const CConfig *config);

  /*!
   * \brief Set flow surface output field values
   * \param[in] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in,out] config - Definition of the particular problem.
   * \param[in] output - Boolean indicating whether information should be written to screen
   */
  void SetAnalyzeSurface(const CSolver* const* solver, const CGeometry *geometry, CConfig *config, bool output);

  /*!
   * \brief Compute and Set flow species variance output field values
   * \param[in] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in,out] config - Definition of the particular problem.
   * \param[in] Surface_Species_Total - Avg mass fraction of each species on all Marker_Analyze
   * \param[in] Surface_MassFlow_Abs_Total - Massflow on all Marker_Analyze
   * \param[in] Surface_Area_Total - Area of all Marker_Analyze
   */
  void SetAnalyzeSurfaceSpeciesVariance(const CSolver* const*solver, const CGeometry *geometry, CConfig *config,
                                         const su2activematrix& Surface_Species_Total,
                                         const vector<su2double>& Surface_MassFlow_Abs_Total,
                                         const vector<su2double>& Surface_Area_Total);

  /*!
   * \brief Add scalar (turbulence/species) history fields for the Residual RMS (FVMComp, FVMInc, FVMNEMO).
   */
  void AddHistoryOutputFields_ScalarRMS_RES(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) history fields for the max Residual (FVMComp, FVMInc, FVMNEMO).
   */
  void AddHistoryOutputFields_ScalarMAX_RES(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) history fields for the BGS Residual (FVMComp, FVMInc, FVMNEMO).
   */
  void AddHistoryOutputFields_ScalarBGS_RES(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) history fields for the linear solver (FVMComp, FVMInc, FVMNEMO).
   */
  void AddHistoryOutputFieldsScalarLinsol(const CConfig* config);

  /*!
   * \brief Set all scalar (turbulence/species) history field values.
   */
  void LoadHistoryDataScalar(const CConfig* config, const CSolver* const* solver);

  /*!
   * \brief Add scalar (turbulence/species) volume solution fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \note The order of fields in restart files is fixed. Therefore the split-up.
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarSolution(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) volume solution fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \note The order of fields in restart files is fixed. Therefore the split-up.
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarResidual(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) volume primitive fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarPrimitive(const CConfig* config);

  /*!
   * \brief Add scalar (turbulence/species) volume limiter fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarLimiter(const CConfig* config);

  /*!
   * \brief Add flamelet volume source term fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarSource(const CConfig* config);

  /*!
   * \brief Add flamelet volume lookup value fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarLookup(const CConfig* config);

  /*!
   * \brief Add miscellaneous scalar volume fields for a point (FVMComp, FVMInc, FVMNEMO).
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFieldsScalarMisc(const CConfig* config);

  /*!
   * \brief Set all scalar (turbulence/species) volume field values for a point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iPoint - Index of the point.
   */
  void LoadVolumeDataScalar(const CConfig* config, const CSolver* const* solver, const CGeometry* geometry,
                             const unsigned long iPoint);

  /*!
   * \brief Add aerodynamic coefficients as output fields
   * \param[in] config - Definition of the particular problem.
   */
  void AddAerodynamicCoefficients(const CConfig* config);

  /*!
   * \brief  Set the value of the aerodynamic coefficients
   * \param[in] config - Definition of the particular problem.
   * \param[in] flow_solver - The container holding all solution data.
   */
  void SetAerodynamicCoefficients(const CConfig* config, const CSolver* flow_solver);

  /*!
   * \brief Add heat flux coefficients as output fields
   * \param[in] config - Definition of the particular problem.
   */
  void AddHeatCoefficients(const CConfig* config);

  /*!
   * \brief  Set the value of the heat flux coefficients
   * \param[in] config - Definition of the particular problem.
   * \param[in] flow_solver - The container holding all solution data.
   */
  void SetHeatCoefficients(const CConfig* config, const CSolver* flow_solver);

  /*!
   * \brief Add rotating frame coefficients as output fields.
   */
  void AddRotatingFrameCoefficients();

  /*!
   * \brief Set the value of the rotating frame coefficients (CT, CQ and CMerit).
   * \param[in] flow_solver - The container holding all solution data.
   */
  void SetRotatingFrameCoefficients(const CSolver* flow_solver);

  /*!
   * \brief Add CP inverse design output as history fields
   */
  void AddCpInverseDesignOutput();

  /*!
   * \brief Set CP inverse design output field values (and also into the solver).
   * \param[in,out] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCpInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Add nearfield inverse design output as history fields
   */
  void AddNearfieldInverseDesignOutput();

  /*!
   * \brief Set nearfield inverse design output field values (and also into the solver).
   * \param[in,out] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNearfieldInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the custom outputs.
   * \param[in] solver - The container holding all solution data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCustomOutputs(const CSolver* const* solver, const CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Helper for custom outputs, converts variable names to indices and pointers which are then used
   * to evaluate the custom expressions.
   */
  void ConvertVariableSymbolsToIndices(const CPrimitiveIndices<unsigned long>& idx, bool allowSkip,
                                       CustomOutput& output) const;

  /*!
   * \brief Compute value of the Q criteration for vortex idenfitication
   * \param[in] VelocityGradient - Velocity gradients
   * \return Value of the Q criteration at the node
   */
  template<class T>
  su2double GetQCriterion(const T& VelocityGradient) const {

    /*--- Make a 3D copy of the gradient so we do not have worry about nDim ---*/

    su2double Grad_Vel[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

    /*--- Q Criterion Eq 1.2 of HALLER, G. (2005). An objective definition of a vortex.
     Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

    /*--- Components of the strain rate tensor (symmetric) ---*/
    su2double s11 = Grad_Vel[0][0];
    su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
    su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
    su2double s22 = Grad_Vel[1][1];
    su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
    su2double s33 = Grad_Vel[2][2];

    /*--- Components of the spin tensor (skew-symmetric) ---*/
    su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
    su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
    su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

    /*--- Q = ||Omega|| - ||Strain|| ---*/
    su2double Q = 2*(pow(omega12,2) + pow(omega13,2) + pow(omega23,2)) -
      (pow(s11,2) + pow(s22,2) + pow(s33,2) + 2*(pow(s12,2) + pow(s13,2) + pow(s23,2)));

    return Q;
  }

  /*!
   * \brief Returns the axisymmetric factor for a point on a marker.
   */
  template <class GeoNodes>
  inline su2double GetAxiFactor(bool axisymmetric, const GeoNodes& nodes, unsigned long iPoint,
                                unsigned short iMarker) {
    if (!axisymmetric) return 1.0;

    if (nodes.GetCoord(iPoint, 1) > EPS) return 2 * PI_NUMBER * nodes.GetCoord(iPoint, 1);

    for (const auto jPoint : nodes.GetPoints(iPoint)) {
      if (nodes.GetVertex(jPoint, iMarker) >= 0) {
        /*--- Not multiplied by two since we need to half the y coordinate. ---*/
        return PI_NUMBER * nodes.GetCoord(jPoint, 1);
      }
    }
    return 0.0;
  }

  /*!
   * \brief Write information to meta data file
   * \param[in] config - Definition of the particular problem per zone.
   */
  void WriteMetaData(const CConfig *config);

  /*!
   * \brief Write any additional files defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - The container holding all solution data.
   */
  void WriteAdditionalFiles(CConfig *config, CGeometry *geometry, CSolver **solver_container) override;

  /*!
   * \brief Names of the currently active volume output fields/groups requested via RESTART_AVG_FIELDS
   *        (matched the same way as VOLUME_OUTPUT: by group name or by individual field name), i.e.
   *        the fields WriteAveragedFields persists and RestoreAveragedFields expects to find.
   * \param[in] config - Definition of the particular problem per zone.
   * \return Field names, in volumeOutput_List order.
   */
  vector<string> GetRequestedAvgFieldNames(const CConfig *config) const;

  /*!
   * \brief Write all currently active TIME_AVERAGE/BACKSCATTER volume output fields to a small
   *        companion ASCII file (together with the number of samples accumulated so far), so they
   *        can be restored on restart instead of restarting the running averages from scratch.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void WriteAveragedFields(CConfig *config, CGeometry *geometry) override;

  /*!
   * \brief Restore all TIME_AVERAGE/BACKSCATTER volume output fields from the companion file
   *        written by WriteAveragedFields, if present.
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void RestoreAveragedFields(CConfig *config, CGeometry *geometry) override;

  /*!
   * \brief Determines if the the volume output should be written.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Iter - Current iteration index.
   * \param[in] force_writing - boolean that forces writing of volume output
   * \param[in] iFile - index to the file that we need to consider for volume output
   */
  bool WriteVolumeOutput(CConfig *config, unsigned long Iter, bool force_writing, unsigned short iFile) override;
  /*!
   * \brief Write the forces breakdown file
   * \param[in] config - Definition of the particular problem per zone.
   * \param[in] flow_solver - The container holding all solution data.
   */
  void WriteForcesBreakdown(const CConfig *config, const CSolver *flow_solver) const;

  /*!
   * \brief Set the time averaged output fields.
   * \param[in] config - Definition of the particular problem per zone.
   */
  void SetTimeAveragedFields(const CConfig *config);

  /*!
   * \brief Load the time averaged output fields.
   * \param iPoint
   * \param node_flow
   * \param node_turb
   * \param config - Definition of the particular problem per zone.
   * \param geometry - Geometrical definition of the problem.
   */
  void LoadTimeAveragedData(unsigned long iPoint, CVariable *node_flow, CVariable *node_turb, const CConfig *config, const CGeometry *geometry);
  /*!
   * \brief Write additional output for fixed CL mode.
   * \param[in] config - Definition of the particular problem per zone.
   */
  void SetFixedCLScreenOutput(const CConfig *config);

  /*!
   * \brief Compute the power of the stochastic forcing (Backscatter Model).
   * \param iPoint - Index of the point.
   * \param config - Definition of the particular problem.
   * \param node_flow - Flow solver solution.
   * \param node_turb - Turbulence-model solver solution.
   * \param geometry - Geometrical definition of the problem.
   * \return Power of the stochastic forcing.
   */
  inline su2double GetPowerStochForcing(unsigned long iPoint, const CConfig *config, const CVariable *node_flow, const CVariable *node_turb, const CGeometry *geometry) {

    const su2double mag = config->GetSBSParam().SBS_Cmag;
    const su2double threshold = config->GetSBSParam().stochFdThreshold;

    /*--- Compute forcing intensity at point i ---*/

    su2double tkeEstim_i = 0.0;
    const su2double lesSensor_i = node_flow->GetLES_Mode(iPoint);

    if (IsHybridRANSLES_SST(config->GetKind_HybridRANSLES())) {
      su2double tke_i = (config->GetSBSParam().useMeanTurb) ? node_turb->GetMeanTurbKinEnergy(iPoint) : node_turb->GetSolution(iPoint, 0);
      tkeEstim_i = (lesSensor_i > threshold) ? tke_i : 0.0;
    } else {
      const su2double lengthscale_i = config->GetConst_DES() * node_turb->GetDES_FilterWidth(iPoint);
      const su2double yoshizawaConst = 0.0066;
      tkeEstim_i = (lesSensor_i > threshold) ? yoshizawaConst * pow(lengthscale_i*node_flow->GetStrainMag(iPoint), 2) : 0.0;
    }

    su2double stochVec_i[3] = {0.0};
    unsigned short startVar = IsHybridRANSLES_SST(config->GetKind_HybridRANSLES()) ? 2 : 1;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      if (config->GetSBSParam().stochSourceType == LANGEVIN)
        stochVec_i[iDim] = node_turb->GetSolution(iPoint, iDim+startVar);
      else if (config->GetSBSParam().stochSourceType == ORNSTEIN_UHLENBECK)
        stochVec_i[iDim] = node_turb->GetOU_Process(iPoint, iDim);
      else
        stochVec_i[iDim] = node_turb->GetLangevinSourceTerms(iPoint, iDim);
      stochVec_i[iDim] *= tkeEstim_i * mag;
    }

    /*--- Evaluate the curl of the stochastic vector ---*/

    su2double curlStochVec[3] = {0.0};

    for (unsigned short iNode = 0; iNode < geometry->nodes->GetnPoint(iPoint); iNode++) {
      auto jPoint = geometry->nodes->GetPoint(iPoint, iNode);
      auto iEdge = geometry->nodes->GetEdge(iPoint, iNode);
      su2double sign = (geometry->edges->GetNode(iEdge, 0) == iPoint) ? 1.0 : -1.0;
      auto* normal = geometry->edges->GetNormal(iEdge);

      /*--- Compute forcing intensity at point j ---*/

      su2double tkeEstim_j = 0.0;
      const su2double lesSensor_j = node_flow->GetLES_Mode(jPoint);

      if (IsHybridRANSLES_SST(config->GetKind_HybridRANSLES())) {
        su2double tke_j = (config->GetSBSParam().useMeanTurb) ? node_turb->GetMeanTurbKinEnergy(jPoint) : node_turb->GetSolution(jPoint, 0);
        tkeEstim_j = (lesSensor_j > threshold) ? tke_j : 0.0;
      } else {
        const su2double lengthscale_j = config->GetConst_DES() * node_turb->GetDES_FilterWidth(jPoint);
        const su2double yoshizawaConst = 0.0066;
        tkeEstim_j = (lesSensor_j > threshold) ? yoshizawaConst * pow(lengthscale_j*node_flow->GetStrainMag(jPoint), 2) : 0.0;
      }

      su2double stochVec_j[3] = {0.0};
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        if (config->GetSBSParam().stochSourceType == LANGEVIN)
          stochVec_j[iDim] = node_turb->GetSolution(jPoint, iDim+startVar);
        else if (config->GetSBSParam().stochSourceType == ORNSTEIN_UHLENBECK)
          stochVec_j[iDim] = node_turb->GetOU_Process(jPoint, iDim);
        else
          stochVec_j[iDim] = node_turb->GetLangevinSourceTerms(jPoint, iDim);
        stochVec_j[iDim] *= tkeEstim_j * mag;
      }

      /*--- Compute fluxes ---*/

      su2double Mean_stochVec[3] = {0.0};
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        Mean_stochVec[iDim] = 0.5 * (stochVec_i[iDim] + stochVec_j[iDim]);
      }

      curlStochVec[0] += sign*(normal[1]*Mean_stochVec[2] - normal[2]*Mean_stochVec[1]);
      curlStochVec[1] += sign*(normal[2]*Mean_stochVec[0] - normal[0]*Mean_stochVec[2]);
      curlStochVec[2] += sign*(normal[0]*Mean_stochVec[1] - normal[1]*Mean_stochVec[0]);
    }

    /*--- Compute the power of the stochastic forcing ---*/

    su2double forcingPower = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      forcingPower += node_flow->GetVelocity(iPoint, iDim) * curlStochVec[iDim];
    }
    forcingPower /= (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint));
    
    return forcingPower;
  }

  /*!
   * \brief Compute the stochastic energy backscatter.
   * \param iPoint - Index of the point.
   * \param config - Definition of the particular problem.
   * \param node_flow - Flow solver solution.
   * \param node_turb - Turbulence-model solver solution.
   * \param geometry - Geometrical definition of the problem.
   * \return Stochastic energy backscatter.
   */
  inline su2double GetEnergyBackscatter(unsigned long iPoint, const CConfig *config, const CVariable *node_flow, const CVariable *node_turb, const CGeometry *geometry) {

    const su2double mag = config->GetSBSParam().SBS_Cmag;
    const su2double threshold = config->GetSBSParam().stochFdThreshold;

    /*--- Compute forcing intensity at point i ---*/

    su2double tkeEstim_i = 0.0;
    const su2double lesSensor_i = node_flow->GetLES_Mode(iPoint);

    if (IsHybridRANSLES_SST(config->GetKind_HybridRANSLES())) {
      su2double tke_i = (config->GetSBSParam().useMeanTurb) ? node_turb->GetMeanTurbKinEnergy(iPoint) : node_turb->GetSolution(iPoint, 0);
      tkeEstim_i = (lesSensor_i > threshold) ? tke_i : 0.0;
    } else {
      const su2double lengthscale_i = config->GetConst_DES() * node_turb->GetDES_FilterWidth(iPoint);
      const su2double yoshizawaConst = 0.0066;
      tkeEstim_i = (lesSensor_i > threshold) ? yoshizawaConst * pow(lengthscale_i*node_flow->GetStrainMag(iPoint), 2) : 0.0;
    }

    su2double stochVec_i[3] = {0.0};
    unsigned short startVar = IsHybridRANSLES_SST(config->GetKind_HybridRANSLES()) ? 2 : 1;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      if (config->GetSBSParam().stochSourceType == LANGEVIN)
        stochVec_i[iDim] = node_turb->GetSolution(iPoint, iDim+startVar);
      else if (config->GetSBSParam().stochSourceType == ORNSTEIN_UHLENBECK)
        stochVec_i[iDim] = node_turb->GetOU_Process(iPoint, iDim);
      else
        stochVec_i[iDim] = node_turb->GetLangevinSourceTerms(iPoint, iDim);
      stochVec_i[iDim] *= tkeEstim_i * mag;
    }

    /*--- Evaluate the dot product of the stochastic vector and the vorticity vector ---*/

    const auto vorticity = node_flow->GetVorticity(iPoint);
    su2double energyBackscatter = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      energyBackscatter += stochVec_i[iDim]*vorticity[iDim];

    return energyBackscatter;
  }
};
