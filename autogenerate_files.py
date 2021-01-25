#
#
#  eT - a coupled cluster program
#  Copyright (C) 2016-2021 the authors of eT
#
#  eT is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  eT is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.
#


# Simple script to autogenerate complex files and interfaces in eT

# Usage: call
# $ python3 autogenerate_interfaces.py
# in the main eT directory (where this file is found).

# Written by Alice Balbi and Andreas Skeidsvoll, Sep 2019
# Restructured by Sander Roet, Feb 2020
# Combined in autogenerate_files.py by Alexander C. Paul, Nov 2020
from re import search, sub, findall
from pathlib import Path
from shutil import rmtree


def is_submodule(file_name):
    # List all submodules found in the wavefunction folders
    # This should be defined as the first not-comment line of the F90 file
    file_name = Path(file_name)
    with file_name.open("r") as f:
        for line in f.readlines():
            # Skip comment lines
            if line.startswith("!"):
                continue
            # First non-comment line start with submodule
            elif line.startswith("submodule"):
                return True
            # First non-comment line does not start with submodule
            else:
                return False


def write_license(f):
    f.write(
        "!\n"
        "!\n"
        "!  eT - a coupled cluster program\n"
        "!  Copyright (C) 2016-2021 the authors of eT\n"
        "!\n"
        "!  eT is free software: you can redistribute it and/or modify\n"
        "!  it under the terms of the GNU General Public "
        "License as published by\n"
        "!  the Free Software Foundation, either version 3 "
        "of the License, or\n"
        "!  (at your option) any later version.\n"
        "!\n"
        "!  eT is distributed in the hope that it will be "
        "useful,\n"
        "!  but WITHOUT ANY WARRANTY; without even the "
        "implied warranty of\n"
        "!  MERCHANTABILITY or FITNESS FOR A PARTICULAR "
        "PURPOSE. See the\n"
        "!  GNU General Public License for more details.\n"
        "!\n"
        "!  You should have received a copy of the GNU "
        "General Public License\n"
        "!  along with this program. If not, see "
        "<https://www.gnu.org/licenses/>.\n"
    )


def fix_line_continuation(lines):
    lines_out = []
    lines_without_comments = []
    for line in lines:
        temp_line = line.rstrip(" ")
        temp_line_without_comments = line.partition("!")[0].strip()

        while temp_line_without_comments.endswith("&"):
            try:
                line_addition = next(lines).rstrip(" ")
            except StopIteration:
                break
            temp_line = temp_line + line_addition
            temp_line_without_comments = (
                temp_line_without_comments.strip("&").strip()
                + line_addition.partition("!")[0].strip()
            )
        lines_out.append(temp_line)
        lines_without_comments.append(temp_line_without_comments)

    return lines_out, lines_without_comments


def find_arguments(line):
    #  Make list of arguments in subroutine or function

    argument_list = []

    #  Finds the text all inside brackets in the given line using regex

    all_inside_brackets = findall(r"\((.*?)\)", line)

    for inside_single_brackets in all_inside_brackets:
        arguments_inside_single_brackets = inside_single_brackets.split(",")
        for argument_name in arguments_inside_single_brackets:
            if argument_name.strip() != "":
                argument_list.append(argument_name.strip().lower())

    return argument_list


def write_function_interface(iterator, arguments, outfile):
    #
    for line, line_without_comments in iterator:
        # Write header lines
        if line[:2] == "!!":
            outfile.write(line)

        # Reached end of subroutine
        if search(r"\b" + "end" + r"\b", line_without_comments) and search(
            r"\b" + "subroutine" + r"\b", line_without_comments
        ):
            outfile.write("!\n" + line)
            break

        # Reached end of function

        if search(r"\b" + "end" + r"\b", line_without_comments) and search(
            r"\b" + "function" + r"\b", line_without_comments
        ):
            outfile.write("!\n" + line)
            break

        # Reached implicit none

        if search(r"\b" + "implicit none" + r"\b", line_without_comments):
            outfile.write(line + "!\n")

        # Reached declaration of variable(s)

        if search("::", line_without_comments):
            after_dots = line_without_comments.partition("::")[2].strip()
            variable_list = [
                variable_name.strip() for variable_name in after_dots.split(",")
            ]
            # Look through all declared variables for the declaration
            # of any subroutine or function arguments

            found_argument = False
            for variable_name in variable_list:
                if variable_name.lower() in arguments:
                    found_argument = True
                    break

            if found_argument:
                outfile.write(line)


def autogenerate_interfaces(source_directory):
    # Autogenerate interfaces
    # Written by Alice Balbi and Andreas Skeidsvoll, Sep 2019
    # Restructured by Sander Roet , Feb 2020
    # Converted to use pathlib.Path and f-strings by Sander Roet, Sep 2020
    # This script reads relevant Fortran files in eT, and makes interfaces.

    wavefunction_directory = source_directory / "wavefunctions"

    # Preparation section
    remove_files = set(wavefunction_directory.glob("*/*_interface*.F90"))
    remove_files.update(set(wavefunction_directory.glob("*/*/*_interface*.F90")))

    # Remove all interface files
    for f in remove_files:
        f.unlink()

    # Get all .F90 files in the directory
    f90_files = set(wavefunction_directory.glob("*/*.F90"))
    f90_files.update(set(wavefunction_directory.glob("*/*/*.F90")))

    wave_function_all_submodules = {i for i in f90_files if is_submodule(i)}

    # Interface generation
    for submodule_name in wave_function_all_submodules:
        with submodule_name.open("r") as f:
            # Make a generator
            lines = (i for i in f.readlines())
        lines, lines_without_comments = fix_line_continuation(lines)
        iterator = zip(lines, lines_without_comments)
        # Strip .F90 and add "_interface.F90"
        interface_name = Path(f"{submodule_name.with_suffix('')}_interface.F90")
        interface_file = interface_name.open("w")
        write_license(interface_file)

        for line, line_without_comments in iterator:
            if search(r"\b" + "module" + r"\b", line_without_comments) and (
                search(r"\b" + "subroutine" + r"\b", line_without_comments)
                or search(r"\b" + "function" + r"\b", line_without_comments)
            ):

                interface_file.write("!\n!\n" + line)

                argument_list = find_arguments(line_without_comments)

                # Loop until reaching end of subroutine or function
                write_function_interface(iterator, argument_list, interface_file)

        interface_file.close()


def autogenerate_complex_files(source_directory):
    # Autogenerate complex files
    # Written by Alice Balbi and Andreas Skeidsvoll, Sep 2019
    # Converted to use pathlib.Path and f-strings by Sander Roet, Sep 2020
    # This script reads relevant real Fortran files in eT, and writes complex
    # versions.

    parameters_file_directory = source_directory / "various"
    wavefunction_directory = source_directory / "wavefunctions"
    autogenerated_directory_name = Path("autogenerated_complex_files")

    # Dictionary with wavefunction directory as key and the corresponding class
    # file as value

    wavefunction_classes = {
        "ccs": ["ccs_class.F90"],
        "doubles": ["doubles_class.F90"],
        "ccsd": ["ccsd_class.F90"],
    }

    # Dictionary with wavefunction directory as key and the corresponding
    # complexified file as value

    wavefunction_complexified_submodules = {
        "ccs": [
            "fock_ccs.F90",
            "initialize_destruct_ccs.F90",
            "jacobian_transpose_ccs.F90",
            "multiplier_equation_ccs.F90",
            "omega_ccs.F90",
            "set_get_ccs.F90",
            "zop_ccs.F90",
            "t1_ccs.F90",
        ],
        "doubles": [
            "initialize_destruct_doubles.F90",
            "jacobian_transpose_doubles.F90",
            "omega_doubles.F90",
            "zop_doubles.F90",
        ],
        "ccsd": [
            "initialize_destruct_ccsd.F90",
            "jacobian_transpose_ccsd.F90",
            "multiplier_equation_ccsd.F90",
            "omega_ccsd.F90",
            "set_get_ccsd.F90",
            "zop_ccsd.F90",
        ],
    }

    excluded_procedures = [
        "initialize_transition_densities_ccs",
        "destruct_transition_densities_ccs",
        "initialize_right_excitation_energies_ccs",
        "destruct_right_excitation_energies_ccs",
        "initialize_left_excitation_energies_ccs",
        "destruct_left_excitation_energies_ccs",
        "ao_to_t1_transformation_ccs",
        "ao_to_t1_transformation_ccs_complex",
        "set_excitation_energies_ccs",
    ]

    # Parameter listing section

    parameters_file = (parameters_file_directory / "parameters.F90").open("r")

    parameter_list = []

    # Read until end of file

    with parameters_file as lines:
        while True:

            try:
                line = next(lines)

            except StopIteration:
                break

            line_without_comments = line.partition("!")[0].strip()

            # In case line contains an ampersand, read next line(s)

            while (len(line_without_comments) > 0) and (
                line_without_comments[-1] == "&"
            ):

                try:
                    line_addition = next(lines)

                except StopIteration:
                    break

                line = line + line_addition
                line_without_comments = (
                    line_without_comments.strip("&").strip()
                    + line_addition.partition("!")[0].strip()
                )

            # Reached a real parameter declaration, add the parameter name to the list
            # of parameters

            if (search(r"\b" + "real" + r"\b", line_without_comments)) and (
                search("::", line_without_comments)
            ):
                #
                after_dots = line_without_comments.partition("::")[2].strip()
                parameter_name = after_dots.partition("=")[0].strip()
                parameter_list.append(parameter_name)

    parameters_file.close()

    # Preparation section

    # Remove autogenerated directories

    for directory_name in wavefunction_classes:

        autogenerated_directory = (
            wavefunction_directory / directory_name / autogenerated_directory_name
        )

        if autogenerated_directory.exists():
            rmtree(autogenerated_directory)

    # List all submodules found in the wavefunction folders

    wavefunction_all_submodules = {}

    for directory_name in wavefunction_classes:

        file_names = []

        for file_name in (wavefunction_directory / directory_name).iterdir():

            file_path = wavefunction_directory / directory_name / file_name

            if file_path.is_file() and (
                file_path.name not in wavefunction_classes[directory_name]
            ):
                file_names.append(file_name)

        wavefunction_all_submodules[directory_name] = file_names

    # Create autogenerated directories

    for directory_name in wavefunction_classes:
        (wavefunction_directory / directory_name / autogenerated_directory_name).mkdir()

    # Make a list containing the name of the real type-bound variables that
    # should be changed to complex, by looking for the type-bound variables
    # before 'contains' in the wavefunction class files

    variable_list = []

    for directory_name in wavefunction_classes:

        class_name = wavefunction_classes[directory_name][0]

        class_file = (wavefunction_directory / directory_name / class_name).open("r")

        # Read until 'contains' or end of file

        with class_file as lines:
            while True:

                try:
                    line = next(lines)

                except StopIteration:
                    break

                line_without_comments = line.partition("!")[0].strip()

                # In case line contains an ampersand, read next line(s)

                while (len(line_without_comments) > 0) and (
                    line_without_comments[-1] == "&"
                ):

                    try:
                        line_addition = next(lines)

                    except StopIteration:
                        break

                    line = line + line_addition
                    line_without_comments = (
                        line_without_comments.strip("&").strip()
                        + line_addition.partition("!")[0].strip()
                    )

                # Reached a variable declaration, add the variable name to the list of
                # real type-bound variables

                if (search(r"\b" + "real" + r"\b", line_without_comments)) and (
                    search("::", line_without_comments)
                ):

                    after_dots = line_without_comments.partition("::")[2].strip()
                    varibles_after_dots = after_dots.split(",")

                    for variable in varibles_after_dots:
                        variable_list.append(variable.strip())

                # Reached the end of type-bound variable declarations

                elif "contains" in line_without_comments:
                    break

        class_file.close()

    # Remove duplicates from the list of variables

    variable_list = list(dict.fromkeys(variable_list))

    # Make a list containing the name of the routines that should be changed to
    # complex, by looking inside the wavefunction submodule files

    procedure_list = []

    for directory_name in wavefunction_complexified_submodules:
        for submodule_name in wavefunction_complexified_submodules[directory_name]:

            real_submodule_file = (
                wavefunction_directory / directory_name / submodule_name
            ).open("r")

            # Read until 'contains' or end of file

            with real_submodule_file as lines:
                while True:

                    try:
                        line = next(lines)

                    except StopIteration:
                        break

                    line_without_comments = line.partition("!")[0].strip()

                    # In case line contains an ampersand, read next line(s)

                    while (len(line_without_comments) > 0) and (
                        line_without_comments[-1] == "&"
                    ):

                        try:
                            line_addition = next(lines)

                        except StopIteration:
                            break

                        line = line + line_addition
                        line_without_comments = (
                            line_without_comments.strip("&").strip()
                            + line_addition.partition("!")[0].strip()
                        )

                    # Reached a subroutine, add the subroutine name to the list of routines

                    if search(
                        r"\b" + "module" + r"\b", line_without_comments
                    ) and search(r"\b" + "subroutine" + r"\b", line_without_comments):

                        before_left_bracket = line_without_comments.partition("(")[
                            0
                        ].strip()

                        subroutine_and_directory_name = before_left_bracket.partition(
                            " subroutine "
                        )[2].strip()

                        if subroutine_and_directory_name != "":
                            procedure_list.append(subroutine_and_directory_name)

                        subroutine_name = subroutine_and_directory_name.rstrip(
                            directory_name
                        ).rstrip("_")

                        if subroutine_name != "":
                            procedure_list.append(subroutine_name)

                    # Reached a function, add the function name to the list of routines

                    if search(
                        r"\b" + "module" + r"\b", line_without_comments
                    ) and search(r"\b" + "function" + r"\b", line_without_comments):

                        before_left_bracket = line_without_comments.partition("(")[
                            0
                        ].strip()

                        function_and_directory_name = before_left_bracket.partition(
                            " function "
                        )[2].strip()

                        if function_and_directory_name != "":
                            procedure_list.append(function_and_directory_name)

                        function_name = function_and_directory_name.rstrip(
                            directory_name
                        ).rstrip("_")

                        if function_name != "":
                            procedure_list.append(function_name)

            real_submodule_file.close()

    procedure_list = list(dict.fromkeys(procedure_list))

    # Submodule section

    for directory_name in wavefunction_complexified_submodules:
        for submodule_name in wavefunction_complexified_submodules[directory_name]:

            real_submodule_file = (
                wavefunction_directory / directory_name / submodule_name
            ).open("r")
            complex_submodule_file = (
                wavefunction_directory
                / directory_name
                / autogenerated_directory_name
                / f"{submodule_name.rstrip('.F90')}_complex.F90"
            ).open("w")

            with real_submodule_file as lines:
                while True:

                    try:
                        line = next(lines)

                    except StopIteration:
                        break

                    line_without_comments = line.partition("!")[0].strip()

                    if ("::" in line_without_comments) and (
                        "ddot" in line_without_comments
                    ):
                        #
                        line = next(lines)
                        line_without_comments = line.partition("!")[0].strip()

                    # Skip excluded procedures

                    if search(r"\b" + "module" + r"\b", line_without_comments) and (
                        search(r"\b" + "subroutine" + r"\b", line_without_comments)
                        or search(r"\b" + "function" + r"\b", line_without_comments)
                    ):

                        for procedure_name in excluded_procedures:

                            if search(
                                r"\b" + procedure_name + r"\b", line_without_comments
                            ):

                                while True:

                                    try:
                                        line = next(lines)

                                    except StopIteration:
                                        break

                                    line_without_comments = line.partition("!")[
                                        0
                                    ].strip()

                                    if (
                                        search(
                                            r"\b" + "end" + r"\b", line_without_comments
                                        )
                                        and (
                                            search(
                                                r"\b" + "subroutine" + r"\b",
                                                line_without_comments,
                                            )
                                            or search(
                                                r"\b" + "function" + r"\b",
                                                line_without_comments,
                                            )
                                        )
                                        and search(
                                            r"\b" + procedure_name + r"\b",
                                            line_without_comments,
                                        )
                                    ):

                                        line = next(lines)

                                        while line.strip() == "!":
                                            line = next(lines)

                                        line_without_comments = line.partition("!")[
                                            0
                                        ].strip()
                                        break

                    # Change submodule names

                    if search(r"\b" + "submodule" + r"\b", line_without_comments):
                        line = line.strip() + "_complex\n"

                    # Change variable declarations

                    line = sub(r"\b" + "real" + r"\b", "complex", line)

                    # Change variable names

                    for variable_name in variable_list:
                        line = sub(
                            r"\b" + "wf%" + variable_name + r"\b",
                            "wf%" + variable_name + "_complex",
                            line,
                        )

                    # Change routine names

                    for procedure_name in procedure_list:
                        line = sub(
                            r"\b" + procedure_name + r"\b",
                            procedure_name + "_complex",
                            line,
                        )

                    # Change wf%eri to wf%eri_c

                    line = sub(r"\b" + "wf%eri" + r"\b", "wf%eri_complex", line)

                    # Change BLAS/LAPACK routines

                    line = sub(r"\b" + "dgemm" + r"\b", "zgemm", line)
                    line = sub(r"\b" + "dcopy" + r"\b", "zcopy", line)
                    line = sub(r"\b" + "daxpy" + r"\b", "zaxpy", line)
                    line = sub(r"\b" + "dscal" + r"\b", "zscal", line)
                    line = sub(r"\b" + "dgemv" + r"\b", "zgemv", line)
                    line = sub(r"\b" + "dger" + r"\b", "zgeru", line)

                    # Change to custom zdotu routine,
                    # since the bundled zdotu routine fails on Macs

                    line = sub(r"\b" + "ddot" + r"\b", "our_zdotu", line)

                    # Change energy and dipole moment

                    line = sub(r"\b" + "energy" + r"\b", "energy_complex", line)
                    line = sub(
                        r"\b" + "dipole_moment" + r"\b", "dipole_moment_complex", line
                    )
                    line = sub(
                        r"\b" + "correlation_energy" + r"\b",
                        "correlation_energy_complex",
                        line,
                    )

                    # Change integral and t1 transformation routines

                    line = sub(r"\b" + "get_t1_oei" + r"\b", "get_t1_oei_complex", line)
                    line = sub(
                        r"\b" + "get_g_pqrs_t1" + r"\b", "get_g_pqrs_t1_complex", line
                    )
                    line = sub(
                        r"\b" + "ao_to_t1_transformation" + r"\b",
                        "ao_to_t1_transformation_complex",
                        line,
                    )

                    # Change array utilities

                    line = sub(r"\b" + "zero_array" + r"\b", "zero_array_complex", line)
                    line = sub(
                        r"\b" + "copy_and_scale" + r"\b", "copy_and_scale_complex", line
                    )

                    # Change parameters

                    for parameter_name in parameter_list:
                        if search(
                            r"\b" + parameter_name + r"\b", line_without_comments
                        ):
                            line = sub(
                                r"\b" + parameter_name + r"\b",
                                parameter_name + "_complex",
                                line,
                            )

                    complex_submodule_file.write(line)

            real_submodule_file.close()
            complex_submodule_file.close()


def main(root_dir):
    source_directory = Path(root_dir / "src")
    autogenerate_complex_files(source_directory)
    autogenerate_interfaces(source_directory)


if __name__ == "__main__":
    default_path = Path(__file__).resolve().parent
    main(default_path)
