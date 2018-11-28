import os
import shutil
from morphct.definitions import PROJECT_ROOT


def main():
    cwd = os.getcwd()
    # Create the directories in the cwd
    os.makedirs(os.path.join(cwd, "inputs"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "outputs"), exist_ok=True)
    # Copy the template files from the MorphCT main dir
    shutil.copy(
        os.path.join(PROJECT_ROOT, "templates", "submit.py"),
        os.path.join(cwd, "submit.py"),
    )
    shutil.copy(
        os.path.join(PROJECT_ROOT, "templates", "submit.sh"),
        os.path.join(cwd, "submit.sh"),
    )
    shutil.copy(
        os.path.join(PROJECT_ROOT, "templates", "par.py"), os.path.join(cwd, "par.py")
    )
