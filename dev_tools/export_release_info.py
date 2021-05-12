#!/usr/bin/env python3
# vim:ft=python
#
#  Licensed under the GNU Lesser General Public License 3.0
#
#  Simple script to extract eT release information and
#  write the information to file for automatic release using pipelines
#
#  Written by Eirik F. Kjønstad, 2021
#

import re
import sys
from datetime import date
from pathlib import Path

version_regex = r"#\s+eT\s+v(\d+)\.(\d+)\.(\d+)"
main_directory = Path(__file__).resolve().parents[1]


def main():

    # Open and read changelog.md

    changelog_file_path = main_directory / "changelog.md"

    with changelog_file_path.open("r") as changelog_file:

        changelog_lines = changelog_file.readlines()

        # Find the most recent release tag and notes (the uppermost release
        # in changelog.md) and export information in the form of files
        # (major.txt, minor.txt, patch.txt, description.txt) that will
        # be read when release is created

        for line_num, line in enumerate(changelog_lines):

            match_version = re.match(version_regex, line)

            if match_version:

                major, minor, patch, name = get_and_check_version(match_version)

                print_version(major, minor, patch, name)

                export_version(major, minor, patch)

                notes = get_release_notes(changelog_lines[line_num + 1 :])

                print_release_notes(notes)

                export_description(notes, major, minor, patch, name)

                break


def get_and_check_version(changelog_match_version):

    # Major, minor, patch according to changelog.md

    major = changelog_match_version.group(1)
    minor = changelog_match_version.group(2)
    patch = changelog_match_version.group(3)

    # Check consistency with src/parameters/various.F90,
    # and get the name of the current minor

    major_, minor_, patch_, name = get_version_from_parameters()

    not_consistent = major != major_ or minor != minor_ or patch != patch_

    if not_consistent:
        raise Exception("Versions in changelog.md and in parameters.F90 do not match!")

    return major, minor, patch, name


def get_version_from_parameters():

    parameters_file_path = main_directory / "src" / "various" / "parameters.F90"

    with parameters_file_path.open("r") as parameters_file:

        lines = parameters_file.readlines()

        major = ""
        minor = ""
        patch = ""
        name = ""

        for line in lines:

            major_match = re.match(r".*major_version\s*=\s*(\d+)", line)
            minor_match = re.match(r".*minor_version\s*=\s*(\d+)", line)
            patch_match = re.match(r".*patch_version\s*=\s*(\d+)", line)
            name_match = re.match(r".*version_name\s*=\s*\"(.*)\"", line)

            if major_match:
                major = major_match.group(1)

            if minor_match:
                minor = minor_match.group(1)

            if patch_match:
                patch = patch_match.group(1)

            if name_match:
                name = name_match.group(1)

            # stop reading if everything is extracted
            if major and minor and patch and name:
                break

    return major, minor, patch, name


def print_version(major, minor, patch, name):

    print(f"Preparing release of eT v{major}.{minor}.{patch} {name}")


def get_release_notes(changelog_lines):

    # changelogLines: starts at the line immediately following
    # the header changelog for the current release (### eT vx.x.x)

    # initialize notes as empty string
    notes = ""

    for line in changelog_lines:

        # exit if you find an earlier version
        if re.match(version_regex, line):
            break

        notes += line

    return notes


def print_release_notes(notes):

    print(f"Release notes: \n{notes}")


def export_version(major, minor, patch):

    major_file_path = main_directory / "major.txt"
    minor_file_path = main_directory / "minor.txt"
    patch_file_path = main_directory / "patch.txt"

    with major_file_path.open("w") as major_file:
        major_file.write(major)

    with minor_file_path.open("w") as minor_file:
        minor_file.write(minor)

    with patch_file_path.open("w") as patch_file:
        patch_file.write(patch)


def export_description(notes, major, minor, patch, name):

    description_file_path = main_directory / "description.txt"

    with description_file_path.open("w") as description_file:

        description_file.write(f"**eT {name}**\\")

        description_file.write(f"Version {major}.{minor}.{patch}\\")

        textual_date = date.today().strftime("%B %d, %Y")  # e.g., "May 08, 2021"

        description_file.write(f"{textual_date}")

        description_file.write(f"{notes}")


if __name__ == "__main__":

    main()
