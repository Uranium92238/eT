#!/usr/bin/env python3
"""

   eT - a coupled cluster program
   Copyright (C) 2016-2021 the authors of eT

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
import shutil
from argparse import ArgumentParser
from pathlib import Path


if sys.version < "3.6":
    print("requires python version >= 3.6")
    sys.exit(1)

default_path = Path.cwd()


def input_parser():
    parser = ArgumentParser()

    parser.add_argument(
        "path",
        help="restart file directory",
        default=default_path,
        action="store",
        nargs="?",
    )

    args = parser.parse_args()

    return args


def path_resolver(args):
    """
    Resolve path to the restart directory.
    """
    translate_path = Path(args.path).resolve()

    print(f"Path provided for converting restart files: '{translate_path}'")

    # Check if path exists
    path_exists = translate_path.is_dir()
    if not path_exists:
        print("Path does not exist")

    return translate_path


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


def get_reference_dimensionalities(file_, endian):
    """
    Read the scf_restart_file to get the number of amplitudes.
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

        bytes_ = f.read(4)
        ee = []
        for i in range(n_states):
            ee.append(f.read(8))

    return ee, n_states


def make_new_gs_file(file_path, energy, n1, n2, endian):
    """
    Make new ground state file.
    n1, n2 are the number of singles and doubles amplitudes.
    Endian specifies the byte order (big or little).
    """
    n1_b = n1.to_bytes(8, byteorder=endian)

    n_records = 2
    if n2:
        n2_b = n2.to_bytes(8, byteorder=endian)
    else:
        n_records = 1

    file_size = file_path.stat().st_size

    # Each record has 4 bytes at the beginning and end
    # Have 1 record if we only have singles, 2 else
    if file_size != 8 * (n1 + n2) + 4 * 2 * n_records:
        raise ValueError(
            f"\nExpected file size {8*(n1+n2+n_records)} found {file_size}.\n"
            f"Number of amplitudes does not match in file {file_path.name}."
        )

    temp = Path("temp")
    with temp.open("w+b") as temp_file:
        temp_file.write(energy)
        temp_file.write(n1_b)

        with file_path.open("rb") as data:
            t1 = read_record(data, endian)
            temp_file.write(t1)

            if n2 > 0:
                temp_file.write(n2_b)
                t2 = read_record(data, endian)
                temp_file.write(t2)

    shutil.move(temp, file_path)


def make_new_es_file(file_path, energy, n1, n2, endian):
    """
    Make new excited state file.
    n1, n2 are the number of singles and doubles amplitudes.
    Endian specifies the byte order (big or little).
    """
    n1_b = n1.to_bytes(8, byteorder=endian)
    if n2:
        n2_b = n2.to_bytes(8, byteorder=endian)

    file_size = file_path.stat().st_size
    # es files only have one record
    if file_size != 8 * (n1 + n2) + 4 * 2:
        raise ValueError(
            f"\nExpected file size {8*(n1+n2+1)} found {file_size}.\n"
            f"Number of amplitudes does not match in file {file_path.name}."
        )

    temp = Path("temp")
    with temp.open("w+b") as temp_file:
        temp_file.write(energy)
        temp_file.write(n1_b)

        with file_path.open("rb") as data:
            # skip record specifier
            bytes_ = data.read(4)
            t1 = data.read(n1 * 8)
            temp_file.write(t1)

            if n2 > 0:
                temp_file.write(n2_b)
                t2 = data.read(8 * n2)
                temp_file.write(t2)

    shutil.move(temp, file_path)


def update_gs_files(dir_, n1, n2, n_gs, endian):
    """
    Update t and tbar file, consider that CC2 ground state files
    have only singles amplitudes.
    Old Format: t1, t2
    New format: 8byte, n_1, singles, n_2, doubles
    """
    file_ = Path(dir_ / "t")

    if file_.exists():
        if n_gs == n1:  # for CC2, n_gs = 0 but n_2 != 0
            make_new_gs_file(file_, bytes(8), n1, 0, endian)
        else:
            make_new_gs_file(file_, bytes(8), n1, n2, endian)

    file_ = Path(dir_ / "tbar")

    if file_.exists():
        if n_gs == n1:
            make_new_gs_file(file_, bytes(8), n1, 0, endian)
        else:
            make_new_gs_file(file_, bytes(8), n1, n2, endian)


def update_es_files(dir_, n1, n2, endian):
    """
    Update excited state files.
    Old Format: full excitation vector in one record
    New format: excitation energy, n_1, singles, n_2, doubles
    """
    ee, n_states = read_excitation_energies(dir_, endian)

    for i in range(1, n_states + 1):
        file_ = Path(dir_ / f"r_{i:03}")
        if file_.exists():
            make_new_es_file(file_, ee[i - 1], n1, n2, endian)

        file_ = Path(dir_ / f"l_{i:03}")
        if file_.exists():
            make_new_es_file(file_, ee[i - 1], n1, n2, endian)


def update_reference_files(dir_, n_ao, endian):
    file_orbitals = Path(dir_ / "orbital_coefficients")
    file_energies = Path(dir_ / "orbital_energies")

    if file_orbitals.exists() and file_energies.exists():
        make_new_orbital_file(file_orbitals, file_energies, n_ao, endian)


def make_new_orbital_file(file_path_orbitals, file_path_energies, n_ao, endian):

    n_ao_b = n_ao.to_bytes(8, byteorder=endian)

    n_records = 1

    file_size = file_path_orbitals.stat().st_size

    if (file_size - 4 * 2 * n_records) > 8 * n_ao ** 2:
        n_records = 2

    # Each record has 4 bytes at the beginning and end
    n_mo = int((file_size - 4 * 2 * n_records) / (n_ao * 8 * n_records))
    n_mo_b = n_mo.to_bytes(8, byteorder=endian)

    temp = Path("temp")
    with temp.open("w+b") as temp_file:
        temp_file.write(n_ao_b)
        temp_file.write(n_mo_b)

        if n_records == 1:

            with file_path_energies.open("rb") as data:
                energies = read_record(data, endian)
                temp_file.write(energies)

            with file_path_orbitals.open("rb") as data:
                orbitals = read_record(data, endian)
                temp_file.write(orbitals)
        else:

            with file_path_energies.open("rb") as data_1:
                energies_a = read_record(data_1, endian)
                energies_b = read_record(data_1, endian)

            with file_path_orbitals.open("rb") as data_2:
                orbitals_a = read_record(data_2, endian)
                orbitals_b = read_record(data_2, endian)

            temp_file.write(energies_a)
            temp_file.write(orbitals_a)
            temp_file.write(energies_b)
            temp_file.write(orbitals_b)

    shutil.move(temp, file_path_orbitals)


def main(argv):
    """
    Read path to restart directory.
    Read scf_restart_file and update orbital files.
    Read cc_restart_file and update ground state files
    Read excitation_energies_file and update excited state files if present.
    """
    args = input_parser()
    translate_path = path_resolver(args)

    # Is the smallest memory address the least significant byte (little endain)
    # or is the smallest address the most significat byte (big endian)?
    endian = sys.byteorder

    restart_file = Path(translate_path / "scf_restart_file")

    if restart_file.exists():
        print("Updating Hartree-Fock files.")
        n_ao = get_reference_dimensionalities(restart_file, endian)
        update_reference_files(translate_path, n_ao, endian)
    else:
        print(
            f"scf_restart_file not found in {translate_path}.\n",
            "Did not update Hartree-Fock restart files.",
        )

    restart_file = Path(translate_path / "cc_restart_file")

    if restart_file.exists():
        n1, n2, n_gs = get_dimensionalities(restart_file, endian)

        if Path(translate_path / "t").exists():
            print("Updating ground state files.")
            update_gs_files(translate_path, n1, n2, n_gs, endian)

        if Path(translate_path / "excitation_energies").exists():
            print("Updating excited state files.")
            update_es_files(translate_path, n1, n2, endian)

    else:
        print(
            f"cc_restart_file not found in {translate_path}.\n",
            "Did not update coupled cluster restart files.",
        )


if __name__ == "__main__":
    main(sys.argv)
