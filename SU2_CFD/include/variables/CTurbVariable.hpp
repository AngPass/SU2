/*!
 * \file CTurbVariable.hpp
 * \brief Base class for defining the variables of the turbulence model.
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

#include "CScalarVariable.hpp"

/*!
 * \class CTurbVariable
 * \brief Base class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbVariable : public CScalarVariable {
protected:
  VectorType muT; /*!< \brief Eddy viscosity. */

public:
  static constexpr size_t MAXNVAR = 5;
  VectorType turb_index;
  VectorType intermittency;         /*!< \brief Value of the intermittency for the trans. model. */

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbVariable() override = default;

  /*!
   * \brief Get the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \return the value of the eddy viscosity.
   */
  inline su2double GetmuT(unsigned long iPoint) const final { return muT(iPoint); }

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \param[in] val_muT - Value of the eddy viscosity.
   */
  inline void SetmuT(unsigned long iPoint, su2double val_muT) final { muT(iPoint) = val_muT; }

  /*!
    * \brief Set the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \param[in] val_turb_index - Value of the turbulence index.
   */
  inline void SetTurbIndex(unsigned long iPoint, su2double val_turb_index) final { turb_index(iPoint) = val_turb_index; }

  /*!
   * \brief Get the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the turbulence index.
   */
  inline su2double GetTurbIndex(unsigned long iPoint) const final { return turb_index(iPoint); }

  /*!
   * \brief Get the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the transition model.
   */
  inline su2double GetIntermittency(unsigned long iPoint) const final { return intermittency(iPoint); }

  /*!
   * \brief Set the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \param[in] val_intermittency - New value of the intermittency.
   */
  inline void SetIntermittency(unsigned long iPoint, su2double val_intermittency) final { intermittency(iPoint) = val_intermittency; }

  /*!
   * \brief Set the Diffusion Coefficients of TKE and omega equations.
   * \param[in] iPoint - Point index.
   * \param[in] val_DC_kw - diffusion coefficient value
   */

  /*!
   * \brief Register eddy viscosity (muT) as Input or Output of an AD recording.
   * \param[in] input - Boolean whether In- or Output should be registered.
   */
  void RegisterEddyViscosity(bool input);

    /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   */
  inline virtual su2double GetLangevinSourceTermsOld(unsigned long iPoint, unsigned short iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_stochSource_old - Old value of source term in Langevin equations.
   */
  inline virtual void SetLangevinSourceTermsOld(unsigned long iPoint, unsigned short iDim, su2double val_stochSource_old) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   */
  inline virtual su2double GetSmoothingMatrixCoeff(unsigned long iPoint, unsigned short iNode) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \param[in] val_smoothmat - IJ-coefficient of the system matrix for Laplacian smoothing of the stochastic source term.
   */
  inline virtual void SetSmoothingMatrixCoeff(unsigned long iPoint, unsigned short iNode, su2double val_smoothmat) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \param[in] iDim - Dimension index.
   */
  inline virtual su2double GetSmoothingBetaVec(unsigned long iPoint, unsigned short iNode, unsigned short iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iNode - Neighbor index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_betavec - Non-orthogonal correction coefficient for the gradient-corrected
   *                           Laplacian smoothing of the stochastic source term.
   */
  inline virtual void SetSmoothingBetaVec(unsigned long iPoint, unsigned short iNode, unsigned short iDim, su2double val_betavec) {}

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   */
  inline virtual su2double GetLangevinSourceGrad(unsigned long iPoint, unsigned short iDim) const { return 0.0; }

  /*!
   * \brief A virtual member.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Dimension index.
   * \param[in] val_grad - Component of the Green-Gauss gradient of the stochastic source term
   *                        currently being smoothed.
   */
  inline virtual void SetLangevinSourceGrad(unsigned long iPoint, unsigned short iDim, su2double val_grad) {}

  /*!
   * \brief Mark the points where the Stochastic Backscatter Model is active.
   */
  inline virtual void SetSBSInBox(unsigned long iPoint, su2double val_sbsInBox) {}

  /*!
   * \brief Get the the points where the Stochastic Backscatter Model is active.
   * \return One if the Stochastic Backscatter Model is active.
   */
  inline virtual su2double GetSBSInBox(unsigned long iPoint) const { return 0.0; }
};
