from pathlib import Path
from tempfile import TemporaryFile
from datetime import datetime
import re
import sys

current_year = datetime.today().year


def allowed(f):
    # List of unreadable suffixes
    unreadable = {".zip", ".swp", ".pyc", ".jpg", ".bin", ".ipynb"}
    excluded_dirs = {"build", ".git", "submodule"}
    if (
        f.is_file()
        and f.suffix not in unreadable
        and not any(str(p) in excluded_dirs for p in f.parents)
    ):
        return True
    return False


def replace_copyright(of):
    # Assume file has been read, so seek to 0
    of.seek(0)
    lines = of.readlines()
    with TemporaryFile("w+") as tf:
        for line in lines:
            ll = line.lower()
            if "copyright" in ll and "et" in ll:
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
                print(f"Updating file: {f}")
                replace_copyright(of)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main()
