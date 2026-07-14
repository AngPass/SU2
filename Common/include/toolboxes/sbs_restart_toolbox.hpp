/*!
 * \file sbs_restart_toolbox.hpp
 * \brief Utilities to persist the per-point, time-averaged ("TIME_AVERAGE"/"BACKSCATTER" volume
 *        output group) fields used by the Stochastic Backscatter Model across restarts, in a small
 *        companion binary file separate from the main restart file, together with the number of
 *        samples already accumulated so the running average can continue coherently.
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

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include "../basic_types/datatype_structure.hpp"
#include "../geometry/CGeometry.hpp"
#include "../parallelization/mpi_structure.hpp"

namespace SBSRestartToolbox {
/// \addtogroup SBSRestartToolbox
/// @{

/*!
 * \brief Binary layout written by WriteMeanFields and read by ReadMeanFields:
 *        uint64 nSamples;
 *        uint64 nFields;
 *        repeated nFields times: uint64 nameLength; char name[nameLength];
 *        repeated once per point (any order): uint64 globalPointIndex; passivedouble value[nFields];
 *        Values are always written/read as passivedouble (i.e. the AD value stripped of derivative
 *        information, same convention used by the native SU2 restart file writers), regardless of
 *        the datatype SU2 was built with. The file is not portable across endianness/precision. ---*/

/*!
 * \brief Write a set of named per-point fields to a small companion binary file, one record per
 *        owned point, keyed by global point index. Ranks write their own points in turn
 *        (round-robin), mirroring the pattern used by CSU2FileWriter for the native restart format.
 * \param[in] filename - Full path of the file to write (overwritten if it exists).
 * \param[in] geometry - Geometry, used to look up global point indices.
 * \param[in] nPointDomain - Number of points owned by this rank.
 * \param[in] fieldNames - Names of the fields being written, in row order.
 * \param[in] nSamples - Number of samples already accumulated in the averages being written.
 * \param[in] fillRow - Callback that fills a row (sized fieldNames.size()) for a given local point.
 */
inline void WriteMeanFields(const std::string& filename, const CGeometry* geometry, unsigned long nPointDomain,
                            const std::vector<std::string>& fieldNames, unsigned long nSamples,
                            const std::function<void(unsigned long, su2double*)>& fillRow) {
  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();
  const auto nFields = fieldNames.size();

  if (rank == 0) {
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    const uint64_t nSamples64 = nSamples;
    file.write(reinterpret_cast<const char*>(&nSamples64), sizeof(nSamples64));
    const uint64_t nFields64 = nFields;
    file.write(reinterpret_cast<const char*>(&nFields64), sizeof(nFields64));
    for (const auto& name : fieldNames) {
      const uint64_t nameLength = name.size();
      file.write(reinterpret_cast<const char*>(&nameLength), sizeof(nameLength));
      file.write(name.data(), nameLength);
    }
  }

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      std::ofstream file(filename, std::ios::binary | std::ios::app);
      std::vector<su2double> row(nFields, 0.0);
      std::vector<passivedouble> passiveRow(nFields);
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        fillRow(iPoint, row.data());
        const uint64_t globalIndex = geometry->nodes->GetGlobalIndex(iPoint);
        file.write(reinterpret_cast<const char*>(&globalIndex), sizeof(globalIndex));
        for (size_t iVar = 0; iVar < nFields; iVar++) passiveRow[iVar] = SU2_TYPE::GetValue(row[iVar]);
        file.write(reinterpret_cast<const char*>(passiveRow.data()), nFields * sizeof(passivedouble));
      }
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

/*!
 * \brief Read a companion mean-fields file written by WriteMeanFields, if present, and apply the
 *        records owned by this rank via the given callback. Every rank reads the whole file
 *        independently and keeps only the points it owns; this is simple and robust, at the cost of
 *        some redundant I/O, acceptable since this happens once at startup for a small file.
 * \param[in] filename - Full path of the file to read.
 * \param[in] geometry - Geometry, used to map global point indices back to local ones.
 * \param[out] outNSamples - Number of samples already accumulated, read from the file header.
 * \param[in] applyRow - Callback invoked with (local point index, field name, value) for owned
 *                        points; the caller decides what to do with each named field (e.g. look up
 *                        its offset in the volume output data and store the value there).
 * \return true if the file was found and read, false if it does not exist.
 */
inline bool ReadMeanFields(const std::string& filename, const CGeometry* geometry, unsigned long& outNSamples,
                           const std::function<void(unsigned long, const std::string&, su2double)>& applyRow) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) return false;

  uint64_t nSamples64 = 0;
  file.read(reinterpret_cast<char*>(&nSamples64), sizeof(nSamples64));
  outNSamples = nSamples64;

  uint64_t nFields64 = 0;
  file.read(reinterpret_cast<char*>(&nFields64), sizeof(nFields64));
  std::vector<std::string> fieldNames(nFields64);
  for (auto& name : fieldNames) {
    uint64_t nameLength = 0;
    file.read(reinterpret_cast<char*>(&nameLength), sizeof(nameLength));
    name.resize(nameLength);
    if (nameLength > 0) file.read(&name[0], nameLength);
  }
  if (fieldNames.empty()) return true;

  const auto nFields = fieldNames.size();
  std::vector<passivedouble> row(nFields);
  uint64_t globalIndex = 0;
  while (file.read(reinterpret_cast<char*>(&globalIndex), sizeof(globalIndex))) {
    if (!file.read(reinterpret_cast<char*>(row.data()), nFields * sizeof(passivedouble))) break;
    const long iPointLocal = geometry->GetGlobal_to_Local_Point(globalIndex);
    if (iPointLocal > -1) {
      for (size_t iVar = 0; iVar < nFields; iVar++)
        applyRow(static_cast<unsigned long>(iPointLocal), fieldNames[iVar], row[iVar]);
    }
  }
  return true;
}

/// @}
}  // namespace SBSRestartToolbox
