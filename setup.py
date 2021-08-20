#!/usr/bin/env python3
# vim:ft=python
#
# Simple setup script to create build directory and run CMake
# Inspired by Dalton setup script by Radovan Bast and Jonas Juselius
# (c) Rolf H. Myhre <rolf.h.myhre at ntnu.no>
# licensed under the GNU Lesser General Public License
# Converted to use pathlib.Path and f-strings by
# Sander Roet <sander.roet at ntnu.no>, Sep 2020
# Added pFUnit options and installation,
# Eirik F. Kjønstad <eirik.kjonstad at ntnu.no>, May 2021

from pathlib import Path

root_dir = Path(__file__).resolve().parent
dev_tool_dir = root_dir / "dev_tools"
default_path = root_dir / "build"
default_pfunit_path = root_dir / "submodules" / "pFUnit" / "build" / "installed"

import sys

sys.path.append(str(dev_tool_dir))

import subprocess
import shutil

from autogenerate_files import autogenerate
from argparse import ArgumentParser, HelpFormatter


class SmartFormatter(HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith("R|"):
            return text[2:].splitlines() + [""]
        # this is the RawTextHelpFormatter._split_lines
        return HelpFormatter._split_lines(self, text, width) + [""]


if sys.version < "3.6":
    print("requires python version >= 3.6")
    sys.exit(1)


def input_parser():
    parser = ArgumentParser(formatter_class=SmartFormatter)

    parser.add_argument(
        "build_dir",
        help="Build directory",
        default=default_path,
        action="store",
        nargs="?",
    )
    parser.add_argument(
        "-clean",
        "--clean-dir",
        help=("Delete build dir and create new (if it exists) " "before running CMake"),
        action="store_true",
    )
    parser.add_argument(
        "-cleanC",
        "--clean-CMake",
        help=("Delete CMakeCache.txt (if it exists) " "before running CMake"),
        action="store_true",
    )

    intgroup = parser.add_mutually_exclusive_group()
    intgroup.add_argument(
        "--int64",
        help="enable 64-bit integers (default)",
        action="store_true",
        default=True,
    )
    intgroup.add_argument(
        "--int32", help="enable 32-bit integers", action="store_true", default=False
    )

    parser.add_argument("-ll", "--libint-lib", help="Path to libint library")

    parser.add_argument("-FC", "--Fortran-compiler", help="Fortran compiler")
    parser.add_argument("-CXX", "--CXX-compiler", help="C++ compiler")
    parser.add_argument("-CC", "--C-compiler", help="C compiler")

    parser.add_argument(
        "-rcheck",
        "--runtime-checks",
        help=("Add runtime checks flag to Fortran compiler " "(debug)"),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-fbatch",
        "--forced-batching",
        help=("Performs batching everywhere " "regardless of memory (debug)"),
        action="store_true",
        default=False,
    )

    ompgroup = parser.add_mutually_exclusive_group()
    ompgroup.add_argument(
        "--omp",
        help="enable OpenMP parallelization (default)",
        action="store_true",
        default=True,
    )
    ompgroup.add_argument(
        "--no-omp",
        help="disable OpenMP parallelization",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--blas-type",
        help="Specify BLAS type to search for.",
        choices=["ATLAS", "MKL", "ACML", "SYSTEM_NATIVE", "ESSL", "OPENBLAS"],
        action="store",
        default=None,
    )

    parser.add_argument(
        "-F-flags",
        "--extra-F-flags",
        help="R|Additional compiler dependent Fortran flags.\n"
        "Note: specify flags as an explicit string starting\n"
        "with whitespace to avoid confusion with setup flags.\n"
        'Example: " -fcoarray=lib"',
    )
    parser.add_argument(
        "-CXX-flags",
        "--extra-CXX-flags",
        help="R|Additional compiler dependent C++ flags.\n"
        "Note: specify flags as an explicit string starting\n"
        "with whitespace to avoid confusion with setup flags.\n"
        'Example: " -Ofast"',
    )
    parser.add_argument(
        "-cmake-flags",
        "--extra-cmake-flags",
        help="R|Additional CMake flags.\n"
        "Note: specify flags as an explicit string starting\n"
        "with whitespace to avoid confusion with setup flags.\n"
        'Example: " -GNinja"',
    )
    parser.add_argument(
        "-linker-flags",
        "--extra-linker-flags",
        help="R|Additional compiler dependent linker flags.\n"
        "Note: specify flags as an explicit string starting\n"
        "with whitespace to avoid confusion with setup flags.\n"
        'Example: " -lcaf_mpi"',
    )
    parser.add_argument(
        "--enable-pfunit",
        help="R|Enable unit testing (pFUnit). Assumes default\n"
        "directory unless custom directory is specified \n"
        "with --pfunit-dir.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--install-pfunit",
        help="R|Auto-install unit testing package (pFUnit).\n"
        "This also enables pFUnit.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--pfunit-dir",
        help="R|Path to pFUnit installation.",
        default=default_pfunit_path,
    )

    # Link to optional libraries
    parser.add_argument(
        "--pcm", help="Enable and link to PCMSolver", action="store_true", default=False
    )

    args = parser.parse_args()

    if args.int32:
        args.int64 = False
    if args.no_omp:
        args.omp = False

    return args


def build_maker(args):
    # Save path to build directory
    build_path = Path(args.build_dir).resolve()

    # Check if user is trying to do something stupid
    eT_dirs = {
        root_dir / "ao_basis",
        root_dir / "cmake",
        root_dir / "documentation",
        root_dir / "src",
        root_dir / "submodules",
        root_dir / "tests",
    }
    if build_path in eT_dirs:
        print(f"'{build_path}' is an eT directory")
        sys.exit(1)

    # Check if it exists
    build_exists = build_path.is_dir()

    # Path to CMakeCache.txt
    Cachename = build_path / "CMakeCache.txt"

    if build_exists:  # Does build already exist
        if args.clean_dir:  # Delete if clean-dir
            shutil.rmtree(build_path)

        # Delete CMakeCache.txt if clean_CMake and it exists
        elif args.clean_CMake:
            if Cachename.exists():
                Cachename.unlink()

    # Make directory if it does not exist
    build_path.mkdir(parents=True, exist_ok=True)

    return build_path


def cmake_translator(arg):
    if arg:
        return "ON"
    else:
        return "OFF"


def construct_command(args):
    command = "cmake"

    if args.install_pfunit:
        args.enable_pfunit = True

    flags = (
        f" -DENABLE_64BIT_INTEGERS={cmake_translator(args.int64)}"
        f" -DENABLE_OMP={cmake_translator(args.omp)}"
        f" -DENABLE_THREADED_MKL={cmake_translator(args.omp)}"
        f" -DENABLE_RUNTIME_CHECKS={cmake_translator(args.runtime_checks)}"
        f" -DENABLE_FORCED_BATCHING={cmake_translator(args.forced_batching)}"
        f" -DENABLE_PCMSOLVER={cmake_translator(args.pcm)}"
        f" -DENABLE_PFUNIT={cmake_translator(args.enable_pfunit)}"
    )
    command += flags

    if args.libint_lib:
        command += f" -DLIBINT2_ROOT={args.libint_lib}"

    if args.Fortran_compiler:
        command += f" -DCMAKE_Fortran_COMPILER={args.Fortran_compiler}"
    if args.CXX_compiler:
        command += f" -DCMAKE_CXX_COMPILER={args.CXX_compiler}"
    if args.C_compiler:
        command += f" -DCMAKE_C_COMPILER={args.C_compiler}"

    if args.extra_F_flags:
        command += f' -DEXTRA_Fortran_FLAGS="{args.extra_F_flags.strip()}"'
    if args.extra_CXX_flags:
        command += f' -DEXTRA_CXX_FLAGS="{args.extra_CXX_flags.strip()}"'
    if args.extra_linker_flags:
        command += f" -DEXTRA_LINKER_FLAGS=" f'"{args.extra_linker_flags.strip()}"'
    if args.extra_cmake_flags:
        command += f" {args.extra_cmake_flags.strip()}"

    if args.blas_type:
        command += f" -DBLAS_TYPE={args.blas_type}"

    if args.enable_pfunit:

        if args.install_pfunit:
            auto_install_pfunit(args)

        command += f" -DCMAKE_PREFIX_PATH={args.pfunit_dir}"

    command += f" {root_dir}"
    return command


def run_CMake(command, build_path):

    output = run_subprocess(command, build_path)
    print(output)

    file_path = build_path / "setup_cmake_output"

    with file_path.open("w") as f:
        f.write(command + "\n\n")
        f.write(output)


def eT_setup(argv):

    args = input_parser()

    autogenerate(root_dir)

    build_path = build_maker(args)
    cmake_command = construct_command(args)
    run_CMake(cmake_command, build_path)

    if args.enable_pfunit:
        copy_pfunit_tests(args, build_path)


def auto_install_pfunit(args):

    pfunit_build = root_dir / "submodules" / "pFUnit" / "build"
    pfunit_install = pfunit_build / "installed"

    print("Will attempt to install pFUnit (path: " + str(pfunit_install) + ")")

    cmake_command = get_pfunit_cmake_command(args, pfunit_build)

    pfunit_build.mkdir(exist_ok=True)

    print(f"Configuring pFUnit: {cmake_command}")

    output = run_subprocess(cmake_command, pfunit_build)
    print(output)

    print("Building and installing pFUnit. This may take a few minutes...")

    output = run_subprocess("make", pfunit_build)
    print(output)

    output = run_subprocess("make install", pfunit_build)
    print(output)


def get_pfunit_cmake_command(args, pfunit_build):

    cmake_command = "cmake -DSKIP_MPI=TRUE "

    skip_omp = "YES"
    if args.omp:
        skip_omp = "NO"

    cmake_command += f"-DSKIP_OPENMP={skip_omp} "

    if args.Fortran_compiler:
        cmake_command += f"-DCMAKE_Fortran_COMPILER={args.Fortran_compiler} "
    if args.C_compiler:
        cmake_command += f"-DCMAKE_C_COMPILER={args.C_compiler} "

    pfunit_src = pfunit_build.resolve().parent

    cmake_command = cmake_command + str(pfunit_src)

    return cmake_command


def copy_pfunit_tests(args, build_path):

    fromDirectory = root_dir / "unit_tests"
    toDirectory = build_path

    shutil.copytree(str(fromDirectory), str(toDirectory), dirs_exist_ok=True)


def run_subprocess(command, path):

    p = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path
    )

    s = p.stderr.decode("utf8")
    s += p.stdout.decode("utf8")

    return s


if __name__ == "__main__":
    eT_setup(sys.argv)
