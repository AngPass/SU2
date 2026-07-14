/*!
 * \file CTurbSSTVariable.cpp
 * \brief Definition of the solution fields.
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


#include "../../include/variables/CTurbSSTVariable.hpp"


CTurbSSTVariable::CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint, unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  sstParsedOptions = config->GetSSTParsedOptions();

  bool backscatter = config->GetSBSParam().StochasticBackscatter;
  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = kine;
    Solution(iPoint,1) = omega;
    if (backscatter && config->GetSBSParam().stochSourceType == LANGEVIN) {
      for (unsigned short iVar = 2; iVar < nVar; iVar++) {
        Solution(iPoint, iVar) = 0.0;
      }
    }
  }

  Solution_Old = Solution;

  sigma_om2 = constants[3];
  beta_star = constants[6];
  prod_lim_const = constants[10];

  F1.resize(nPoint) = su2double(1.0);
  F2.resize(nPoint) = su2double(0.0);
  CDkw.resize(nPoint) = su2double(0.0);

  muT.resize(nPoint) = mut;

  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
    DES_LengthScale.resize(nPoint) = su2double(0.0);
    DES_FilterWidth.resize(nPoint) = su2double(0.0);
    MeanTurbKE.resize(nPoint) = su2double(0.0);
    lesMode.resize(nPoint) = su2double(0.0);
    Vortex_Tilting.resize(nPoint) = su2double(0.0);
    if (backscatter) {
      stochSource.resize(nPoint, nDim) = su2double(0.0);
      stochSourceOld.resize(nPoint, nDim) = su2double(0.0);
      besselIntegral.resize(nPoint) = su2double(0.0);
      OU_Process.resize(nPoint, nDim) = su2double(0.0);
      sbsInBox.resize(nPoint) = su2double(1.0);
      smoothMatrix.resize(nPoint, MAXNNEIGHBORS) = su2double(0.0);
      smoothBetaVec.resize(nPoint, MAXNNEIGHBORS*3) = su2double(0.0);
      langevinSourceGrad.resize(nPoint, 3) = su2double(0.0);
      smoothDiag.resize(nPoint) = su2double(1.0);
      smoothPhat.resize(nPoint, 3) = su2double(0.0);
      smoothShat.resize(nPoint, 3) = su2double(0.0);
    }
  }
}

void CTurbSSTVariable::SetVortex_Tilting(unsigned long iPoint, CMatrixView<const su2double> PrimGrad_Flow,
                                         const su2double* Vorticity, su2double LaminarViscosity) {

  su2double Strain[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}}, Omega, StrainDotVort[3], numVecVort[3];
  su2double numerator, trace0, trace1, denominator;

  AD::StartPreacc();
  AD::SetPreaccIn(PrimGrad_Flow, nDim+1, nDim);
  AD::SetPreaccIn(Vorticity, 3);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(muT(iPoint));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(LaminarViscosity);

  Strain[0][0] = PrimGrad_Flow[1][0];
  Strain[1][0] = 0.5*(PrimGrad_Flow[2][0] + PrimGrad_Flow[1][1]);
  Strain[0][1] = 0.5*(PrimGrad_Flow[1][1] + PrimGrad_Flow[2][0]);
  Strain[1][1] = PrimGrad_Flow[2][1];
  if (nDim == 3){
    Strain[0][2] = 0.5*(PrimGrad_Flow[3][0] + PrimGrad_Flow[1][2]);
    Strain[1][2] = 0.5*(PrimGrad_Flow[3][1] + PrimGrad_Flow[2][2]);
    Strain[2][0] = 0.5*(PrimGrad_Flow[1][2] + PrimGrad_Flow[3][0]);
    Strain[2][1] = 0.5*(PrimGrad_Flow[2][2] + PrimGrad_Flow[3][1]);
    Strain[2][2] = PrimGrad_Flow[3][2];
  }

  Omega = sqrt(Vorticity[0]*Vorticity[0] + Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);

  StrainDotVort[0] = Strain[0][0]*Vorticity[0]+Strain[0][1]*Vorticity[1]+Strain[0][2]*Vorticity[2];
  StrainDotVort[1] = Strain[1][0]*Vorticity[0]+Strain[1][1]*Vorticity[1]+Strain[1][2]*Vorticity[2];
  StrainDotVort[2] = Strain[2][0]*Vorticity[0]+Strain[2][1]*Vorticity[1]+Strain[2][2]*Vorticity[2];

  numVecVort[0] = StrainDotVort[1]*Vorticity[2] - StrainDotVort[2]*Vorticity[1];
  numVecVort[1] = StrainDotVort[2]*Vorticity[0] - StrainDotVort[0]*Vorticity[2];
  numVecVort[2] = StrainDotVort[0]*Vorticity[1] - StrainDotVort[1]*Vorticity[0];

  numerator = sqrt(6.0) * sqrt(numVecVort[0]*numVecVort[0] + numVecVort[1]*numVecVort[1] + numVecVort[2]*numVecVort[2]);
  trace0 = 3.0*(pow(Strain[0][0],2.0) + pow(Strain[1][1],2.0) + pow(Strain[2][2],2.0));
  trace1 = pow(Strain[0][0] + Strain[1][1] + Strain[2][2],2.0);
  denominator = pow(Omega, 2.0) * sqrt(trace0-trace1);

  Vortex_Tilting(iPoint) = (numerator/denominator) * max(1.0,0.2*LaminarViscosity/muT(iPoint));

  AD::SetPreaccOut(Vortex_Tilting(iPoint));
  AD::EndPreacc();
}

void CTurbSSTVariable::SetBlendingFunc(unsigned long iPoint, su2double val_viscosity,
                                       su2double val_dist, su2double val_density, TURB_TRANS_MODEL trans_model) {
  su2double arg2, arg2A, arg2B, arg1;

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Solution[iPoint], nVar);
  AD::SetPreaccIn(Gradient[iPoint], nVar, nDim);

  /*--- Cross diffusion ---*/

  CDkw(iPoint) = 0.0;
  for (unsigned long iDim = 0; iDim < nDim; iDim++)
    CDkw(iPoint) += Gradient(iPoint,0,iDim)*Gradient(iPoint,1,iDim);
  CDkw(iPoint) *= 2.0*val_density*sigma_om2/Solution(iPoint,1);
  CDkw(iPoint) = max(CDkw(iPoint), pow(10.0, -prod_lim_const));

  /*--- F1 ---*/

  arg2A = sqrt(Solution(iPoint,0))/(beta_star*Solution(iPoint,1)*val_dist+EPS*EPS);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution(iPoint,1)+EPS*EPS);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution(iPoint,0) / (CDkw(iPoint)*val_dist*val_dist+EPS*EPS));
  F1(iPoint) = tanh(pow(arg1, 4.0));

  /*--- F2 ---*/

  arg2 = max(2.0*arg2A, arg2B);
  F2(iPoint) = tanh(pow(arg2, 2.0));

  /*--- LM model for F1 ---*/
  if (trans_model == TURB_TRANS_MODEL::LM) {
    su2double Ry = val_density*val_dist*sqrt(Solution(iPoint,0))/val_viscosity;
    su2double F3 = exp(-pow(Ry/120.0, 8.0));
    F1(iPoint) = max(F1(iPoint), F3);
  }

  AD::SetPreaccOut(F1(iPoint)); AD::SetPreaccOut(F2(iPoint)); AD::SetPreaccOut(CDkw(iPoint));
  AD::EndPreacc();

}
