#!/usr/bin/env python

## \file seed_restart_averages.py
#  \brief Seed the WRT_RESTART_AVERAGES companion file (RESTART_FILENAME_average_XXXXX.dat) from a
#         steady RANS solution or an existing mean-velocity field (restart/solution file or .vtu), so
#         an unsteady run (DDES/URANS) can start FILTER_STRESSES with a physically-reasonable mean
#         flow instead of accumulating it from scratch.
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
Builds a WRT_RESTART_AVERAGES companion file directly from an existing mean-velocity field (a steady
RANS restart/solution file, or a .vtu volume output snapshot), instead of accumulating it by running
an unsteady simulation first.

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

Which fields can be derived, and from which source
----------------------------------------------------
--solution accepts three kinds of input, auto-detected by content/extension:
  - A native SU2 ASCII restart/solution file (.csv, or .dat without the binary magic number): derives
    MEAN_VELOCITY-X/Y/Z from "Velocity_x/y/z" (incompressible) or "Momentum_x/y/z" / "Density"
    (compressible), and MEAN_TURB_KIN_ENERGY from "Turb_Kin_Energy" (SST only).
  - A native SU2 binary restart/solution file (.dat, with the binary magic number): same fields, same
    derivation, just parsed from the binary layout instead.
  - A SU2 volume output file (.vtu, ParaView XML unstructured grid), e.g. a snapshot written mid-run
    with VOLUME_OUTPUT including TIME_AVERAGE fields: reads MEAN_VELOCITY-X/Y/Z directly from the
    "MeanVelocity" point-data vector (SU2 bundles same-named "_x"/"_y"/"_z" scalar fields into one
    vector array for ParaView, see CParaviewXMLFileWriter::WriteData in the C++ source -- so
    MEAN_VELOCITY-X/Y/Z appears there as a single 3-component array named "MeanVelocity"), and
    MEAN_TURB_KIN_ENERGY from the scalar "MeanTurbulentKineticEnergy". Parsed by a small,
    dependency-free reader (no VTK/meshio required) written specifically against SU2's own .vtu output
    (uncompressed, "appended"/"raw" encoding -- the only variant SU2 ever writes; a .vtu re-saved or
    re-encoded by another tool, e.g. to base64 or ASCII inline data, is not supported). Use
    --velocity-field/--tke-field if your .vtu uses different array names (e.g. a plain steady RANS
    .vtu with just "Velocity", or a field renamed by another tool). Reading a .vtu also loses
    precision: SU2 writes volume output fields as float32, vs. float64 in a native restart file --
    fine for seeding a mean that will keep evolving, less so if you need bit-for-bit reproducibility.
  IMPORTANT for all three: the point order in the source file must match the mesh's global point
  index order (point i in the file = global point index i), which holds for a full-volume file written
  by SU2 itself, but NOT for a decimated/surface-only export or one reordered by another tool.
These are exactly the RESTART_AVG_FIELDS SU2 defaults to when WRT_RESTART_AVERAGES=YES (or, with
WRT_RESTART_AVERAGES=NO, when RESTART_AVERAGE=YES and FILTER_STRESSES is active -- the "frozen mean"
case) and RESTART_AVG_FIELDS is left empty. If you have customized RESTART_AVG_FIELDS to include
additional fields, this script will NOT seed those (SU2 will not error, but each missing field
silently restarts from zero while being weighted as if it already had --samples prior samples -- see
the warning SU2 prints on restore, and CFlowOutput::RestoreAveragedFields in the C++ source).

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

Examples
--------
    python seed_restart_averages.py --solution solution_flow.csv --output solution_flow_average_00000.dat
    python seed_restart_averages.py --solution flow_snapshot.vtu --output solution_flow_average_00000.dat
"""

import argparse
import struct
import sys
import xml.etree.ElementTree as ET

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
        help="Path to the source file: a converged steady RANS restart/solution file (.csv/.dat ASCII, "
        "or .dat binary), or a SU2 volume output file (.vtu) containing a mean/instantaneous velocity "
        "field. See the 'Which fields can be derived' section in this script's module docstring.",
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
    parser.add_argument(
        "--velocity-field",
        default=None,
        help="Base name of the velocity point-data array to use instead of the built-in candidates "
        "(MeanVelocity, Velocity). For a .vtu, this is the bundled vector array's name (e.g. "
        "'Velocity' if the array is named 'Velocity' rather than split into Velocity_x/y/z); for a "
        "restart/solution file, the scalar columns '<name>_x'/'<name>_y'/'<name>_z' are used.",
    )
    parser.add_argument(
        "--tke-field",
        default=None,
        help="Name of the turbulent kinetic energy point-data array to use instead of the built-in "
        "candidates (MeanTurbulentKineticEnergy, Turb_Kin_Energy).",
    )
    args = parser.parse_args()

    if args.samples < 1:
        parser.error("--samples must be at least 1.")

    header, data = read_su2_solution(args.solution)
    fields, values = derive_mean_fields(header, data, args.velocity_field, args.tke_field)

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
    """Read a SU2 native restart/solution file (ASCII .csv, or binary .dat) or a SU2 volume output
    file (.vtu), auto-detected. Returns (header, data) where header is a list of column names
    (including "PointID") and data is a 2D numpy array of shape (nPoints, nColumns) in the same
    column order as header."""

    if path.lower().endswith(".vtu"):
        return _read_vtu(path)

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


_VTK_TYPE_TO_NUMPY = {
    "Float32": "f4", "Float64": "f8",
    "Int8": "i1", "Int16": "i2", "Int32": "i4", "Int64": "i8",
    "UInt8": "u1", "UInt16": "u2", "UInt32": "u4", "UInt64": "u8",
}


def _read_vtu(path):
    """Parse a .vtu written by SU2's own writer (CParaviewXMLFileWriter::WriteData), by hand, without
    any third-party VTK library: this format is a small, fixed, fully-documented subset of the VTK XML
    UnstructuredGrid spec (single Piece, "appended"/"raw" encoding, no compression -- SU2 never writes
    any other variant), so a dependency-free reader tailored to exactly this is both simpler and more
    robust than pulling in a general-purpose VTK reader for the one thing we need: the PointData
    arrays. (A third-party library is more general but also has far more format variations to get
    right; a hand-rolled parser scoped to SU2's own fixed writer output has none of that surface.)"""

    with open(path, "rb") as f:
        raw = f.read()

    # The header (plain ASCII XML) and the appended binary blob must be separated before doing any XML
    # parsing, since the blob can contain arbitrary bytes that are not valid XML/UTF-8. Per the VTK XML
    # spec, "appended" data starts with a single "_" right after the AppendedData tag's ">".
    marker = b'encoding="raw">\n_'
    marker_pos = raw.find(marker)
    if marker_pos == -1:
        marker = b'encoding="raw">_'
        marker_pos = raw.find(marker)
    if marker_pos == -1:
        raise ValueError(
            "'{}': could not find a raw/appended AppendedData section. This reader only supports "
            "SU2's own .vtu format (uncompressed, appended, raw encoding); a .vtu re-saved by another "
            "tool (e.g. ASCII or base64-encoded) is not supported.".format(path)
        )
    head = raw[: marker_pos + len(marker)]
    blob_start = marker_pos + len(marker)
    blob_end = raw.find(b"</AppendedData>", blob_start)
    if blob_end == -1:
        raise ValueError("'{}': could not find the closing </AppendedData> tag.".format(path))
    blob = raw[blob_start:blob_end]

    # Close out the two still-open elements (AppendedData, VTKFile) so the head parses as valid XML;
    # its own content (the lone "_") is irrelevant, we already sliced the real binary data into "blob".
    root = ET.fromstring(head + b"</AppendedData></VTKFile>")

    byte_order = root.attrib.get("byte_order", "LittleEndian")
    endian = "<" if byte_order == "LittleEndian" else ">"
    header_dtype = np.dtype(endian + ("u8" if root.attrib.get("header_type") == "UInt64" else "u4"))
    header_size = header_dtype.itemsize

    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError("'{}': no <Piece> element found.".format(path))
    nPoints = int(piece.attrib["NumberOfPoints"])

    point_data = piece.find("PointData")
    if point_data is None:
        raise ValueError("'{}': no <PointData> found (nothing to read).".format(path))

    header = ["PointID"]
    columns = [np.arange(nPoints, dtype=np.float64)]
    suffixes = ["x", "y", "z"]

    for data_array in point_data.findall("DataArray"):
        name = data_array.attrib["Name"]
        nComp = int(data_array.attrib.get("NumberOfComponents", "1"))
        offset = int(data_array.attrib["offset"])
        vtk_type = data_array.attrib.get("type", "Float32")
        if vtk_type not in _VTK_TYPE_TO_NUMPY:
            raise ValueError("'{}': unsupported DataArray type '{}' for '{}'.".format(path, vtk_type, name))
        dtype = np.dtype(endian + _VTK_TYPE_TO_NUMPY[vtk_type])

        # Per-array layout (CParaviewXMLFileWriter::WriteDataArray): a header-sized byte count,
        # immediately followed by that many bytes of raw, tightly-packed data, both starting at the
        # array's declared "offset" from the first byte of the blob (i.e. right after the "_" marker).
        nbytes = int(np.frombuffer(blob, dtype=header_dtype, count=1, offset=offset)[0])
        values = np.frombuffer(blob, dtype=dtype, count=nbytes // dtype.itemsize, offset=offset + header_size)
        values = values.astype(np.float64)
        if nPoints * nComp != values.size:
            raise ValueError(
                "'{}': array '{}' has {} value(s), expected {} ({} points x {} component(s)).".format(
                    path, name, values.size, nPoints * nComp, nPoints, nComp
                )
            )

        if nComp > 1:
            values = values.reshape(nPoints, nComp)
            for i in range(nComp):
                header.append("{}_{}".format(name, suffixes[i]))
                columns.append(values[:, i])
        else:
            header.append(name)
            columns.append(values)

    return header, np.column_stack(columns)


# Candidate point-data array names to look for, in priority order, for each derived field. A restart
# file lists the velocity components as separate scalar columns ("<name>_x"/"_y"/"_z"); a .vtu read via
# _read_vtu has already been unbundled into the same shape, so the same candidates work for both.
VELOCITY_FIELD_CANDIDATES = ["MeanVelocity", "Velocity"]
TKE_FIELD_CANDIDATES = ["MeanTurbulentKineticEnergy", "Turb_Kin_Energy"]


def derive_mean_fields(header, data, velocity_field=None, tke_field=None):
    """Derive the MEAN_VELOCITY-X/Y/Z and MEAN_TURB_KIN_ENERGY columns from the columns/point-data
    actually present, trying (in order) the explicit --velocity-field/--tke-field override if given,
    then the built-in candidate names, then (for velocity only) the compressible-restart fallback of
    deriving it from Momentum_i/Density. Returns (field_names, field_value_arrays), only including
    fields that could actually be derived."""

    col = {name: idx for idx, name in enumerate(header)}
    fields, values = [], []

    for name in ([tke_field] if tke_field else TKE_FIELD_CANDIDATES):
        if name in col:
            fields.append("MEAN_TURB_KIN_ENERGY")
            values.append(data[:, col[name]])
            break

    velocity = None
    vel_bases = [velocity_field] if velocity_field else VELOCITY_FIELD_CANDIDATES
    for base in vel_bases:
        vx, vy, vz = base + "_x", base + "_y", base + "_z"
        if vx in col and vy in col:
            velocity = [data[:, col[vx]], data[:, col[vy]]]
            if vz in col:
                velocity.append(data[:, col[vz]])
            break

    if velocity is None and "Momentum_x" in col and "Density" in col:
        # Compressible native restart file: velocity_i = momentum_i / density.
        rho = data[:, col["Density"]]
        velocity = [data[:, col["Momentum_x"]] / rho, data[:, col["Momentum_y"]] / rho]
        if "Momentum_z" in col:
            velocity.append(data[:, col["Momentum_z"]] / rho)

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
