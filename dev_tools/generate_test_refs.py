#
#
#  eT - a coupled cluster program
#  Copyright (C) 2016-2022 the authors of eT
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
"""
Copies all .out files not ending in .timing.out from the subdirectories 
under eT/build/tests to corresponding subdirectories in eT/tests.
Obviously, you have to make sure all tests are passing before running this script.
It will modify all reference files, so it will be obvious in a diff.
Written by Rolf H. Myhre, Mar 2021
"""

from pathlib import Path
from argparse import ArgumentParser


def copy_test(to_path, from_path):
    """
    Copy all files with .out suffix except files ending in .timing.out
    from from_path to to_path and remove information sensitive to national security.
    """

    outpaths = [
        f
        for f in from_path.iterdir()
        if f.suffix == ".out" and f.suffixes[0] != ".timing"
    ]

    for outpath in outpaths:

        with outpath.open("r") as outfile:
            lines = outfile.readlines()

        copy_path = to_path / "result" / outpath.name

        with copy_path.open("w") as copy_file:
            for line in lines:
                copy_file.write(line)


def main(root_path):

    # Set up paths to eT/tests and eT/build/tests
    test_path = root_path / "tests"
    build_path = root_path / "build" / "tests"

    # List all directories in eT/tests/
    testlist = [test.name for test in test_path.iterdir() if test.is_dir()]

    for test in testlist:
        copy_test(test_path / test, build_path / test)


if __name__ == "__main__":

    # Take eT root directory from input or use calling directory as default
    parser = ArgumentParser()
    parser.add_argument("eT_dir", help="eT base directory", nargs="?", metavar="eT dir")
    args = parser.parse_args()
    if args.eT_dir:
        eT_path = Path(args.eT_dir).resolve()
    else:
        eT_path = Path(__file__).resolve().parents[1]

    main(eT_path)
