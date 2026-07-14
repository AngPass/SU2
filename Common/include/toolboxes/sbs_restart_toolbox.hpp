/*!
 * \file sbs_restart_toolbox.hpp
 * \brief Utilities to persist the per-point, time-averaged ("TIME_AVERAGE"/"BACKSCATTER" volume
 *        output group) fields used by the Stochastic Backscatter Model across restarts, in a small
 *        companion ASCII file separate from the main restart file, together with the number of
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

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include "../geometry/CGeometry.hpp"
#include "../parallelization/mpi_structure.hpp"

namespace SBSRestartToolbox {
/// \addtogroup SBSRestartToolbox
/// @{

/*!
 * \brief Write a set of named per-point fields to a small companion ASCII file, one row per owned
 *        point, keyed by global point index. Ranks write their own points in turn (round-robin),
 *        mirroring the pattern used by CSU2FileWriter for the native restart format. The first line
 *        stores the number of samples already accumulated in the averages being written.
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
    std::ofstream file(filename);
    file << "# N_SAMPLES " << nSamples << "\n";
    file << "\"PointID\"";
    for (const auto& name : fieldNames) file << ",\"" << name << "\"";
    file << "\n";
  }

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      std::ofstream file(filename, std::ios::app);
      file.precision(15);
      std::vector<su2double> row(nFields, 0.0);
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        fillRow(iPoint, row.data());
        file << geometry->nodes->GetGlobalIndex(iPoint);
        for (const auto& value : row) file << ", " << std::scientific << value;
        file << "\n";
      }
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

/*!
 * \brief Read a companion mean-fields file written by WriteMeanFields, if present, and apply the
 *        rows owned by this rank via the given callback. Every rank reads the whole file
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
  std::ifstream file(filename);
  if (!file.is_open()) return false;

  std::string line;

  /*--- First line: "# N_SAMPLES <value>". ---*/
  std::getline(file, line);
  {
    std::istringstream iss(line);
    std::string hash, tag;
    iss >> hash >> tag >> outNSamples;
  }

  /*--- Second line: comma-separated, quoted field names, "PointID" first. ---*/
  std::getline(file, line);
  std::vector<std::string> fieldNames;
  {
    std::istringstream iss(line);
    std::string token;
    while (std::getline(iss, token, ',')) {
      /*--- Strip surrounding quotes. ---*/
      const auto first = token.find('"');
      const auto last = token.rfind('"');
      if (first != std::string::npos && last != std::string::npos && last > first)
        fieldNames.push_back(token.substr(first + 1, last - first - 1));
    }
  }
  if (fieldNames.empty()) return true;
  fieldNames.erase(fieldNames.begin()); /*--- Drop "PointID", the remaining entries are the data fields. ---*/

  while (std::getline(file, line)) {
    if (line.empty()) continue;
    for (char& c : line)
      if (c == ',') c = ' ';
    std::istringstream iss(line);

    unsigned long globalIndex;
    iss >> globalIndex;

    const long iPointLocal = geometry->GetGlobal_to_Local_Point(globalIndex);
    for (const auto& name : fieldNames) {
      su2double value;
      iss >> value;
      if (iPointLocal > -1) applyRow(static_cast<unsigned long>(iPointLocal), name, value);
    }
  }
  return true;
}

/// @}
}  // namespace SBSRestartToolbox
