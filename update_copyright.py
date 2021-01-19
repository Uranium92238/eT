from pathlib import Path
from tempfile import TemporaryFile
from datetime import datetime
import re
import sys

current_year = datetime.today().year


def allowed(f):
    # List of unreadable suffixes
    unreadable = {".zip", ".swp", ".pyc", ".jpg", ".bin", ".ipynb"}
    excluded_dirs = {"build", ".git", "submodules", "docs"}
    if (
        f.is_file()
        and f.suffix not in unreadable
        and not any(str(p) in excluded_dirs for p in f.parents)
    ):
        return True
    return False


def update_autodocs():
    # Special function for updating the copyright of the doc autogeneration
    # script
    f = Path("./documentation/eT_docs.md")
    with f.open("r+") as of:
        # We are looking for the 'date' line here instead of copyright
        replace_copyright(of, requirements={"date"})


def replace_copyright(of, requirements={"copyright", "et"}):
    # Assume file has been read, so seek to 0
    of.seek(0)
    lines = of.readlines()
    with TemporaryFile("w+") as tf:
        for line in lines:
            ll = line.lower()
            if all(req in ll for req in requirements):
                print(f"Updating file: {of.name}")
                # Regex explaination:
                # "20[0-9][0-9]" matches 20**
                # (?! ....) negative lookahead (unmatch a match if .... is found
                # after it).
                # ".*20[0-9][0-9]" [0-inf] characters followed by 20**
                # Effectively this says 20** not followed by 20**, or
                # replace last occurence of the number 20**
                line = re.sub("20[0-9][0-9](?!.*20[0-9][0-9])", str(current_year), line)
            tf.write(line)
        # Now overwrite the original file
        of.seek(0)
        tf.seek(0)
        of.write(tf.read())


def main(rootdir=""):
    print(f"Updating copyright year to {current_year}")
    root = Path("")
    files = [f for f in root.rglob("*") if allowed(f)]
    for f in files:
        with f.open("r+") as of:
            try:
                content = of.read()
            except UnicodeDecodeError:
                continue
            if "copyright" in content.lower():
                replace_copyright(of)
    update_autodocs()


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main()
