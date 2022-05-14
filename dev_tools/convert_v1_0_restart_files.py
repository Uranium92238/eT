#!/usr/bin/env python3
"""

   eT - a coupled cluster program
   Copyright (C) 2016-2022 the authors of eT

   eT is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program. If not, see <https://www.gnu.org/licenses/>.

   Written by Alexander C. Paul and Sarai D. Folkestad, 2020

"""

import sys
from argparse import ArgumentParser
from pathlib import Path


if sys.version_info < (3, 6):
    print("requires python version >= 3.6")
    sys.exit(1)


def input_parser():
    parser = ArgumentParser()

    parser.add_argument(
        "old_restart_dir",
        help="Old restart file directory",
    )

    parser.add_argument(
        "new_restart_dir",
        help="New restart file directory",
        nargs="?",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="Print verbose information",
        action="store_true",
    )

    args = parser.parse_args()

    return args


def path_resolver(args):
    """
    Resolve paths to the restart directories.
    """
    v0_path = Path(args.v0_restart_dir).resolve()
    if args.vx_restart_dir:
        vx_path = Path(args.vx_restart_dir).resolve()
    else:
        vx_path = v0_path / "v_1_x_restart_files"

    if args.verbose:
        print(f"v1.0 restart directory: {str(v0_path)}")
        print(f"v1.x restart directory: {str(vx_path)}")

    # Check that v1.0 path exists and is not the same as v1.x path
    if not v0_path.is_dir():
        raise NotADirectoryError(str(v0_path))
    if v0_path == vx_path:
        raise Exception("v1.0 directory cannot be the same as v1.x directory.")

    return v0_path, vx_path


def read_record(file_, endian):
    """
    Read a single record of the shape:
        4 bytes representing length of the record
        Data
        4 bytes representing length of the record
    Endian specifies the byte order (big or little).
    """
    bytes_ = file_.read(4)
    record_length = int.from_bytes(bytes_, endian)
    record = file_.read(record_length)
    bytes_ = file_.read(4)

    return record


def get_dimensionalities(file_, endian):
    """
    Read the cc_restart_file to get the number of amplitudes.
    Endian specifies the byte order (big or little).
    """
    with file_.open("rb") as f:
        n_o = int.from_bytes(read_record(f, endian), endian)
        n_v = int.from_bytes(read_record(f, endian), endian)
        n_gs = int.from_bytes(read_record(f, endian), endian)
        n_es = int.from_bytes(read_record(f, endian), endian)

        n1 = n_o * n_v
        if n_es == n1:  # CCS
            n2 = 0
        else:
            n2 = n1 * (n1 + 1) // 2

    return n1, n2, n_gs


def get_n_ao(file_, endian):
    """
    Read the scf_restart_file to get the number of atomic orbitals.
    Endian specifies the byte order (big or little).
    """
    with file_.open("rb") as f:
        n_ao = int.from_bytes(read_record(f, endian), endian)

    return n_ao


def read_excitation_energies(dir_, endian):
    """
    Read the excitation_energies file.
    Endian specifies the byte order (big or little).
    """
    file_ = Path(dir_ / "excitation_energies")
    with file_.open("rb") as f:
        n_states = int.from_bytes(read_record(f, endian), endian)

        f.read(4)
        ee = []
        for i in range(n_states):
            ee.append(f.read(8))

    return ee, n_states


def make_new_gs_file(v0_path, vx_path, energy, n1, n2, endian):
    """
    Make new ground state file.
    n1, n2 are the number of singles and doubles amplitudes.
    Endian specifies the byte order (big or little).
    """
    n1_b = n1.to_bytes(8, byteorder=endian)

    if n2 > 0:
        n2_b = n2.to_bytes(8, byteorder=endian)
        n_records = 2
    else:
        n_records = 1

    file_size = v0_path.stat().st_size

    # Each record has 4 bytes at the beginning and end
    # Have 1 record if we only have singles, 2 else
    if file_size != 8 * (n1 + n2) + 4 * 2 * n_records:
        raise ValueError(
            f"\nExpected file size {8*(n1+n2+n_records)} found {file_size}.\n"
            f"Number of amplitudes does not match in file {str(v0_path)}."
        )

    with vx_path.open("wb") as vx_file, v0_path.open("rb") as v0_file:
        vx_file.write(energy)
        vx_file.write(n1_b)

        t1 = read_record(v0_file, endian)
        vx_file.write(t1)

        if n2 > 0:
            vx_file.write(n2_b)
            t2 = read_record(v0_file, endian)
            vx_file.write(t2)


def make_new_es_file(v0_path, vx_path, energy, n1, n2, endian):
    """
    Make new excited state file.
    n1, n2 are the number of singles and doubles amplitudes.
    Endian specifies the byte order (big or little).
    """
    n1_b = n1.to_bytes(8, byteorder=endian)
    if n2:
        n2_b = n2.to_bytes(8, byteorder=endian)

    file_size = v0_path.stat().st_size
    if file_size != 8 * (n1 + n2) + 4 * 2:
        raise ValueError(
            f"\nExpected file size {8*(n1+n2+1)} found {file_size}.\n"
            f"Number of amplitudes does not match in file {str(v0_path)}."
        )

    # es files only have one record
    with vx_path.open("wb") as vx_file, v0_path.open("rb") as v0_file:
        vx_file.write(energy)
        vx_file.write(n1_b)

        # skip record specifier
        v0_file.read(4)
        t1 = v0_file.read(n1 * 8)
        vx_file.write(t1)

        if n2 > 0:
            vx_file.write(n2_b)
            t2 = v0_file.read(8 * n2)
            vx_file.write(t2)


def make_new_cc_gs_files(v0_path, vx_path, n1, n2, n_gs, endian):
    """
    Update t and tbar file, consider that CC2 ground state files
    have only singles amplitudes.
    Old Format: t1, t2
    New format: 8byte, n_1, singles, n_2, doubles
    """
    v0_file = v0_path / "t"
    vx_file = vx_path / "t"

    if v0_file.is_file():
        if n_gs == n1:  # for CC2, n_gs = n1 but n_2 != 0
            make_new_gs_file(v0_file, vx_file, bytes(8), n1, 0, endian)
        else:
            make_new_gs_file(v0_file, vx_file, bytes(8), n1, n2, endian)

    v0_file = v0_path / "tbar"
    vx_file = vx_path / "tbar"

    if v0_file.is_file():
        if n_gs == n1:
            make_new_gs_file(v0_file, vx_file, bytes(8), n1, 0, endian)
        else:
            make_new_gs_file(v0_file, vx_file, bytes(8), n1, n2, endian)


def make_new_cc_es_files(v0_path, vx_path, n1, n2, endian):
    """
    Update excited state files.
    Old Format: full excitation vector in one record
    New format: excitation energy, n_1, singles, n_2, doubles
    """
    ee, n_states = read_excitation_energies(v0_path, endian)

    for i in range(n_states):
        v0_file = v0_path / f"r_{i+1:03}"
        vx_file = vx_path / f"r_{i+1:03}"
        if v0_file.is_file():
            make_new_es_file(v0_file, vx_file, ee[i], n1, n2, endian)

        v0_file = v0_path / f"l_{i+1:03}"
        vx_file = vx_path / f"l_{i+1:03}"
        if v0_file.is_file():
            make_new_es_file(v0_file, vx_file, ee[i], n1, n2, endian)


def make_new_orbital_file(v0_path, vx_path, n_ao, endian):
    """
    Make new orbital file.
    Old format: separate files for energies, coefficients and dimensions
    New format: n_ao, n_mo, energies, coefficients, (energies, coefficients)
    """

    orbitals_file = v0_path / "orbital_coefficients"
    if not orbitals_file.is_file():
        raise FileNotFoundError(str(orbitals_file))

    energies_file = v0_path / "orbital_energies"
    if not energies_file.is_file():
        raise FileNotFoundError(str(energies_file))

    n_ao_b = n_ao.to_bytes(8, byteorder=endian)

    file_size = orbitals_file.stat().st_size

    # Each record has 4 bytes at the beginning and end
    if file_size - 8 <= 8 * n_ao**2:
        n_records = 1
    elif file_size - 16 <= 2 * 8 * n_ao**2:
        n_records = 2
    else:
        raise Exception(
            f"File size of file {str(orbitals_file)} " "does not match n_ao = {n_ao}"
        )

    n_mo = int((file_size - 4 * 2 * n_records) / (n_ao * 8 * n_records))
    n_mo_b = n_mo.to_bytes(8, byteorder=endian)

    vx_path = Path(vx_path / "orbital_coefficients")
    with vx_path.open("wb") as vx_file:
        vx_file.write(n_ao_b)
        vx_file.write(n_mo_b)

        if n_records == 1:

            with energies_file.open("rb") as data:
                energies = read_record(data, endian)
                vx_file.write(energies)

            with orbitals_file.open("rb") as data:
                orbitals = read_record(data, endian)
                vx_file.write(orbitals)
        else:

            with energies_file.open("rb") as data_1, orbitals_file.open("rb") as data_2:
                energies = read_record(data_1, endian)
                vx_file.write(energies)

                orbitals = read_record(data_2, endian)
                vx_file.write(orbitals)

                energies = read_record(data_1, endian)
                vx_file.write(energies)

                orbitals = read_record(data_2, endian)
                vx_file.write(orbitals)


def convert_v1_0_files(v0_path, vx_path, verbose):
    """
    Read scf_restart_file and update orbital files.
    Read cc_restart_file and update ground state files
    Read excitation_energies_file and update excited state files if present.
    """
    # Is the smallest memory address the least significant byte (little endain)
    # or is the smallest address the most significat byte (big endian)?
    endian = sys.byteorder

    v0_path = Path(v0_path).resolve()
    vx_path = Path(vx_path).resolve()

    vx_path.mkdir(parents=True, exist_ok=True)

    restart_file = v0_path / "scf_restart_file"

    if restart_file.is_file():
        if verbose:
            print("Creating new Hartree-Fock files.")
        n_ao = get_n_ao(restart_file, endian)
        make_new_orbital_file(v0_path, vx_path, n_ao, endian)
    elif verbose:
        print(
            f"scf_restart_file not found in {v0_path}.\n",
            "Did not update Hartree-Fock restart files.",
        )

    restart_file = v0_path / "cc_restart_file"
    excitation_file = v0_path / "excitation_energies"

    if restart_file.is_file():
        n1, n2, n_gs = get_dimensionalities(restart_file, endian)

        if verbose:
            print("Creating new CC ground state files.")
        make_new_cc_gs_files(v0_path, vx_path, n1, n2, n_gs, endian)

        if excitation_file.is_file():
            if verbose:
                print("Creating new CC excited state files.")
            make_new_cc_es_files(v0_path, vx_path, n1, n2, endian)

    elif verbose:
        print(
            f"cc_restart_file not found in {v0_path}.\n",
            "Did not update coupled cluster restart files.",
        )


def main():
    """
    Read paths to restart directories and pass on to converter.
    """
    args = input_parser()
    v0_path, vx_path = path_resolver(args)
    convert_v1_0_files(v0_path, vx_path, args.verbose)


if __name__ == "__main__":
    main()
