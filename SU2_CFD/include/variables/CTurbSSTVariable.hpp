/*!
 * \file CTurbSSTVariable.hpp
 * \brief Declaration of the variables of the SST turbulence model.
 * \author F. Palacios, T. Economon
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

#include "CTurbVariable.hpp"

/*!
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSSTVariable final : public CTurbVariable {
protected:
  su2double sigma_om2;
  su2double beta_star;
  su2double prod_lim_const;
  VectorType F1;
  VectorType F2;    /*!< \brief Menter blending function for blending of k-w and k-eps. */
  VectorType CDkw;  /*!< \brief Cross-diffusion. */
  SST_ParsedOptions sstParsedOptions;
  VectorType DES_LengthScale;
  VectorType MeanTurbKE;
  VectorType lesMode;
  MatrixType stochSource;
  MatrixType stochSourceOld;
  MatrixType OU_Process;
  VectorType sbsInBox;
  VectorType besselIntegral;
  MatrixType smoothMatrix;
  MatrixType smoothBetaVec;      /*!< \brief Non-orthogonal correction coefficients (packed per neighbor as 3 components) for the gradient-corrected Laplacian smoothing. */
  MatrixType langevinSourceGrad; /*!< \brief Green-Gauss gradient of the stochastic source term component currently being smoothed. */
  su2double MAXNNEIGHBORS = 32;
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] kine - Turbulence kinetic energy (k) (initialization value).
   * \param[in] omega - Turbulent variable value (initialization value).
   * \param[in] mut - Eddy viscosity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants - sst model constants
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSSTVariable() override = default;

  /*!
   * \brief Get the DES length scale
   * \param[in] iPoint - Point index.
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(unsigned long iPoint) const override { return DES_LengthScale(iPoint); }

  /*!
   * \brief Set the DES Length Scale.
   * \param[in] iPoint - Point index.
   */
  inline void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) override { DES_LengthScale(iPoint) = val_des_lengthscale; }

  /*!
   * \brief Get the mean turbulent kinetic energy.
   * \param[in] iPoint - Point index.
   * \return Mean turbulent kinetic energy.
   */
  inline su2double GetMeanTurbKinEnergy(unsigned long iPoint) const override { return MeanTurbKE(iPoint); }

  /*!
   * \brief Set the mean turbulent kinetic energy.
   */
  inline void SetMeanTurbKinEnergy(unsigned long iPoint, su2double val_mean_tke) override { MeanTurbKE(iPoint) = val_mean_tke; }

  /*!
   * \brief Get the source terms for the stochastic equations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \return Value of the source term for the stochastic equations.
   */
  inline su2double GetLangevinSourceTerms(unsigned long iPoint, unsigned short iDim) const override { return stochSource(iPoint, iDim); }

  /*!
   * \brief Set the source terms for the stochastic equations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_stochSource - Value of the source term for the stochastic equations.
   */
  inline void SetLangevinSourceTerms(unsigned long iPoint, unsigned short iDim, su2double val_stochSource) override { stochSource(iPoint, iDim) = val_stochSource; }

  /*!
   * \brief Get the old value of the source terms for the stochastic equations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \return Old value of the source terms for the stochastic equations.
   */
  inline su2double GetLangevinSourceTermsOld(unsigned long iPoint, unsigned short iDim) const override { return stochSourceOld(iPoint, iDim); }

  /*!
   * \brief Set the old value of source terms for the stochastic equations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_stochSource_old - Old value of the source term for the stochastic equations.
   */
  inline void SetLangevinSourceTermsOld(unsigned long iPoint, unsigned short iDim, su2double val_stochSource_old) override { stochSourceOld(iPoint, iDim) = val_stochSource_old; }

  /*!
   * \brief Get the IJ-coefficient of the system matrix for Laplacian smoothing of the stochastic source term.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \return IJ-coefficient of the system matrix for Laplacian smoothing of the stochastic source term.
   */
  inline su2double GetSmoothingMatrixCoeff(unsigned long iPoint, unsigned short iNode) const override { return smoothMatrix(iPoint, iNode); }

  /*!
   * \brief Set the IJ-coefficient of the system matrix for Laplacian smoothing of the stochastic source term.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \param[in] val_smoothmat - IJ-coefficient of the system matrix for Laplacian smoothing of the stochastic source term.
   */
  inline void SetSmoothingMatrixCoeff(unsigned long iPoint, unsigned short iNode, su2double val_smoothmat) override { smoothMatrix(iPoint, iNode) = val_smoothmat; }

  /*!
   * \brief Get a component of the non-orthogonal correction coefficient for the gradient-corrected
   *        Laplacian smoothing of the stochastic source term.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \param[in] iDim - Dimension index.
   */
  inline su2double GetSmoothingBetaVec(unsigned long iPoint, unsigned short iNode, unsigned short iDim) const override {
    return smoothBetaVec(iPoint, iNode*3 + iDim);
  }

  /*!
   * \brief Set a component of the non-orthogonal correction coefficient for the gradient-corrected
   *        Laplacian smoothing of the stochastic source term.
   */
  inline void SetSmoothingBetaVec(unsigned long iPoint, unsigned short iNode, unsigned short iDim, su2double val_betavec) override {
    smoothBetaVec(iPoint, iNode*3 + iDim) = val_betavec;
  }

  /*!
   * \brief Get a component of the Green-Gauss gradient of the stochastic source term currently
   *        being smoothed.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   */
  inline su2double GetLangevinSourceGrad(unsigned long iPoint, unsigned short iDim) const override {
    return langevinSourceGrad(iPoint, iDim);
  }

  /*!
   * \brief Set a component of the Green-Gauss gradient of the stochastic source term currently
   *        being smoothed.
   */
  inline void SetLangevinSourceGrad(unsigned long iPoint, unsigned short iDim, su2double val_grad) override {
    langevinSourceGrad(iPoint, iDim) = val_grad;
  }

  /*!
   * \brief Get the value of the stochastic variable resulting from the Ornstein-Uhlenbeck process.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \return Value of the stochastic variable resulting from Ornstein-Uhlenbeck process.
   */
  inline su2double GetOU_Process(unsigned long iPoint, unsigned short iDim) const override { return OU_Process(iPoint, iDim); }

  /*!
   * \brief Set the value of the stochastic variable resulting from the Ornstein-Uhlenbeck process.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_OU - Value of the stochastic variable resulting from Ornstein-Uhlenbeck process.
   */
  inline void SetOU_Process(unsigned long iPoint, unsigned short iDim, su2double val_OU) override { OU_Process(iPoint, iDim) = val_OU; }

  /*!
   * \brief Set the LES sensor.
   */
  inline void SetLES_Mode(unsigned long iPoint, su2double val_les_mode) override { lesMode(iPoint) = val_les_mode; }

  /*!
   * \brief Get the LES sensor.
   * \return Value of the LES sensor.
   */
  inline su2double GetLES_Mode(unsigned long iPoint) const override { return lesMode(iPoint); }

  /*!
   * \brief Mark the points where the Stochastic Backscatter Model is active.
   */
  inline void SetSBSInBox(unsigned long iPoint, su2double val_sbsInBox) override { sbsInBox(iPoint) = val_sbsInBox; }

  /*!
   * \brief Get the the points where the Stochastic Backscatter Model is active.
   * \return One if the Stochastic Backscatter Model is active.
   */
  inline su2double GetSBSInBox(unsigned long iPoint) const override { return sbsInBox(iPoint); }

  /*!
   * \brief Set the integral of the product of three Bessel functions appearing in Laplacian smoothing.
   * \param[in] iPoint - Point index.
   * \param[in] val_integral - Value of the integral.
   */
  inline void SetBesselIntegral(unsigned long iPoint, su2double val_integral) override { besselIntegral(iPoint) = val_integral; }

  /*!
   * \brief Get the the integral of the product of three Bessel functions appearing in Laplacian smoothing.
   * \return Value of the integral.
   */
  inline su2double GetBesselIntegral(unsigned long iPoint) const override { return besselIntegral(iPoint); }

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   */
  void SetBlendingFunc(unsigned long iPoint, su2double val_viscosity, su2double val_dist, su2double val_density, TURB_TRANS_MODEL trans_model) override;

  /*!
   * \brief Get the first blending function.
   */
  inline su2double GetF1blending(unsigned long iPoint) const override { return F1(iPoint); }

  /*!
   * \brief Get the second blending function.
   */
  inline su2double GetF2blending(unsigned long iPoint) const override { return F2(iPoint); }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline su2double GetCrossDiff(unsigned long iPoint) const override { return CDkw(iPoint); }
};
