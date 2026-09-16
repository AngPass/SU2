#!/usr/bin/env python

## \file seed_restart_averages.py
#  \brief Seed the WRT_RESTART_AVERAGES companion file (RESTART_FILENAME_average_XXXXX.dat) from a
#         converged steady RANS solution, so an unsteady run (DDES/URANS) can start FILTER_STRESSES
#         with a physically-reasonable mean flow instead of accumulating it from scratch.
#  \author SU2 Team
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

"""
Builds a WRT_RESTART_AVERAGES companion file directly from a steady RANS restart/solution file,
instead of accumulating it by running an unsteady simulation first.

Background
----------
When WRT_RESTART_AVERAGES=YES, SU2 persists a small companion binary file next to the main restart
file (<SOLUTION_FILENAME>_average_XXXXX.dat, see CFlowOutput::WriteAveragedFields/RestoreAveragedFields
in the C++ source) containing the running time-average of a few fields (by default
MEAN_TURB_KIN_ENERGY and MEAN_VELOCITY-X/Y/Z) plus the number of samples already accumulated. On
restart (RESTART_AVERAGE=YES), this lets the running average continue coherently instead of resetting.

A steady RANS solution IS already an estimate of the mean flow: its velocity field and turbulent
kinetic energy (SST) are, by construction, time-averaged quantities. This script repackages a
converged steady RANS restart/solution file into that same companion-file format, so you can:

    1. Run a normal steady RANS case to convergence.
    2. Run this script on its solution file to produce a *_average_XXXXX.dat file.
    3. Start the unsteady (DDES/URANS) run with RESTART_SOL=YES pointing at the steady solution
       as initial condition, WRT_RESTART_AVERAGES=YES, RESTART_AVERAGE=YES, and
       RESTART_ITER_AVERAGE matching the iteration this script's output is named for (see "Output
       filename" below), so FILTER_STRESSES can be switched on from a mean flow that is already
       close to converged, rather than starting from zero.

Binary format (must match Common/include/toolboxes/sbs_restart_toolbox.hpp exactly):
    uint64 nSamples
    uint64 nFields
    repeated nFields times: uint64 nameLength; char name[nameLength]   (not null-terminated)
    repeated once per point (any order): uint64 globalPointIndex; float64 value[nFields]
All integers/floats are written little-endian, double precision, matching a standard (non-AD) SU2
build. This is the same assumption the C++ writer itself makes ("not portable across endianness or
precision"); it will not match a build using single precision or a different endianness.

Which fields can be derived
----------------------------
From a plain RANS restart file (i.e. the default SOLUTION-group columns, without extra fields added
via VOLUME_OUTPUT), this script derives:
  - MEAN_VELOCITY-X/Y/Z: read directly (incompressible solver, columns "Velocity_x/y/z") or computed
    as Momentum_i / Density (compressible solver, columns "Momentum_x/y/z" and "Density").
  - MEAN_TURB_KIN_ENERGY: read directly from "Turb_Kin_Energy" (SST only; not available for SA, since
    SA has no separate k equation -- the script simply omits this field in that case).
These are exactly the RESTART_AVG_FIELDS SU2 defaults to when WRT_RESTART_AVERAGES=YES and
RESTART_AVG_FIELDS is left empty. If you have customized RESTART_AVG_FIELDS to include additional
fields, this script will NOT seed those (SU2 will not error, but each missing field silently restarts
from zero while being weighted as if it already had --samples prior samples -- see the warning SU2
prints on restore, and CFlowOutput::RestoreAveragedFields in the C++ source).

Output filename
----------------
SU2 looks for <SOLUTION_FILENAME (no extension)>_average<suffix>.dat, where <suffix> is
"_" + the iteration number zero-padded to 5 digits (e.g. "_00000") whenever TIME_DOMAIN=YES, which it
always is for this feature. The iteration used is RESTART_ITER_AVERAGE if set, otherwise RESTART_ITER.
For a fresh unsteady run seeded from a steady solution, RESTART_ITER is typically 0, so the target is
usually <SOLUTION_FILENAME>_average_00000.dat. Point this script's --output at that exact path.

Choosing --samples
-------------------
The seeded file's sample count controls how strongly the running average clings to the steady-RANS
seed once the unsteady run starts (CUMULATIVE averaging blends in each new instantaneous sample at
weight 1/(n+1)):
  - A small value (e.g. 1) lets the true unsteady mean take over quickly -- good if you expect the
    DDES/URANS mean to meaningfully differ from the steady RANS mean (e.g. massively separated flow),
    but does little to shorten the practical accumulation time before FILTER_STRESSES is fully
    "at strength".
  - A large value (e.g. 1000) keeps the average anchored near the RANS seed for a long time -- good if
    you trust the RANS mean and mainly want to skip the numerical transient, but will slow (or mostly
    suppress) any real drift of the mean away from the RANS solution.
There is no universally correct value; pick it based on how much you trust the steady RANS mean for
your case. The default below (100) is a moderate middle ground, not a validated recommendation.

Example
-------
    python seed_restart_averages.py --solution solution_flow.csv --output solution_flow_average_00000.dat
"""

import argparse
import struct
import sys

import numpy as np

DEFAULT_SAMPLES = 100
BINARY_MAGIC = 535532
CGNS_STRING_SIZE = 33


def main():
    parser = argparse.ArgumentParser(
        description="Seed a WRT_RESTART_AVERAGES companion file from a steady RANS restart/solution file."
    )
    parser.add_argument(
        "--solution",
        required=True,
        help="Path to the converged steady RANS restart/solution file (.csv/.dat ASCII, or .dat binary).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the companion average file to (must match what SU2 will look for, see the "
        "'Output filename' section in this script's module docstring).",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=DEFAULT_SAMPLES,
        help="Number of prior samples to record in the seeded file (default: %(default)s). See the "
        "'Choosing --samples' section in this script's module docstring.",
    )
    args = parser.parse_args()

    if args.samples < 1:
        parser.error("--samples must be at least 1.")

    header, data = read_su2_solution(args.solution)
    fields, values = derive_mean_fields(header, data)

    if not fields:
        sys.exit(
            "Error: could not derive any of MEAN_VELOCITY-X/Y/Z or MEAN_TURB_KIN_ENERGY from the "
            "columns in '{}'. Found columns: {}".format(args.solution, ", ".join(header))
        )

    point_ids = data[:, header.index("PointID")].astype(np.uint64)
    write_average_file(args.output, point_ids, fields, values, args.samples)

    print("Wrote {} field(s) for {} point(s), {} prior samples, to '{}'.".format(
        len(fields), len(point_ids), args.samples, args.output))
    print("Fields written: " + ", ".join(fields))


def read_su2_solution(path):
    """Read a SU2 native restart/solution file, ASCII (.csv) or binary (.dat), auto-detected.
    Returns (header, data) where header is a list of column names (including "PointID") and data is a
    2D numpy array of shape (nPoints, nColumns) in the same column order as header."""

    with open(path, "rb") as f:
        first_bytes = f.read(4)

    is_binary = len(first_bytes) == 4 and struct.unpack("<i", first_bytes)[0] == BINARY_MAGIC

    if is_binary:
        return _read_su2_binary(path)
    return _read_su2_ascii(path)


def _read_su2_ascii(path):
    with open(path, "r") as f:
        header_line = f.readline()
    header = [h.strip().strip('"') for h in header_line.strip().split(",")]
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[1] != len(header):
        raise ValueError(
            "'{}': header has {} column(s) but data rows have {}.".format(path, len(header), data.shape[1])
        )
    return header, data


def _read_su2_binary(path):
    with open(path, "rb") as f:
        magic, nVar, nPoint = struct.unpack("<iii", f.read(12))
        f.read(8)  # two reserved/padding ints
        if magic != BINARY_MAGIC:
            raise ValueError("'{}' does not look like a SU2 binary restart file.".format(path))

        header = []
        for _ in range(nVar):
            raw = f.read(CGNS_STRING_SIZE)
            header.append(raw.split(b"\x00", 1)[0].decode("ascii").strip())

        data = np.frombuffer(f.read(8 * nVar * nPoint), dtype="<f8").reshape(nPoint, nVar)

    # Unlike the ASCII writer, the binary format (CSU2BinaryFileWriter::WriteData) has no explicit
    # PointID column: a row's position in the file IS its global point index (rows are written in
    # increasing global-index order). Synthesize the column so both readers return the same shape.
    header = ["PointID"] + header
    data = np.column_stack([np.arange(nPoint, dtype=np.float64), data])
    return header, data


def derive_mean_fields(header, data):
    """Derive the MEAN_VELOCITY-X/Y/Z and MEAN_TURB_KIN_ENERGY columns from a RANS solution's native
    field set. Returns (field_names, field_value_arrays), only including fields that could be derived
    from the columns actually present."""

    col = {name: idx for idx, name in enumerate(header)}
    fields, values = [], []

    if "Turb_Kin_Energy" in col:
        fields.append("MEAN_TURB_KIN_ENERGY")
        values.append(data[:, col["Turb_Kin_Energy"]])

    if "Velocity_x" in col:
        # Incompressible solver: velocity is a primary solution variable.
        velocity = [data[:, col["Velocity_x"]], data[:, col["Velocity_y"]]]
        if "Velocity_z" in col:
            velocity.append(data[:, col["Velocity_z"]])
    elif "Momentum_x" in col and "Density" in col:
        # Compressible solver: velocity_i = momentum_i / density.
        rho = data[:, col["Density"]]
        velocity = [data[:, col["Momentum_x"]] / rho, data[:, col["Momentum_y"]] / rho]
        if "Momentum_z" in col:
            velocity.append(data[:, col["Momentum_z"]] / rho)
    else:
        velocity = None

    if velocity is not None:
        for comp, name in zip(velocity, ["MEAN_VELOCITY-X", "MEAN_VELOCITY-Y", "MEAN_VELOCITY-Z"]):
            fields.append(name)
            values.append(comp)

    return fields, values


def write_average_file(path, point_ids, field_names, field_values, nSamples):
    """Write the companion binary file, matching SBSRestartToolbox::WriteMeanFields exactly
    (Common/include/toolboxes/sbs_restart_toolbox.hpp): little-endian, double precision."""

    nPoints = len(point_ids)
    nFields = len(field_names)

    with open(path, "wb") as f:
        f.write(struct.pack("<Q", nSamples))
        f.write(struct.pack("<Q", nFields))
        for name in field_names:
            name_bytes = name.encode("ascii")
            f.write(struct.pack("<Q", len(name_bytes)))
            f.write(name_bytes)

        record_dtype = np.dtype([("idx", "<u8"), ("vals", "<f8", (nFields,))])
        records = np.zeros(nPoints, dtype=record_dtype)
        records["idx"] = point_ids
        records["vals"] = np.column_stack(field_values)
        f.write(records.tobytes())


if __name__ == "__main__":
    main()
