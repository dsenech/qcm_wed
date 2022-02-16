import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

import os
import subprocess
git_hash = str(subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']))
git_hash = git_hash[2:-3]

fout = open("pyqcm/qcm_git_hash.py", "a")
fout.write("git_hash = '{:s}'\n".format(git_hash))
fout.close() 

setup(
    name="qcm",
    version="1.0",
    description="Quantum cluster methods for the physics of strongly correlated systems",
    author="David Sénéchal",
    license="MIT",
    packages=find_packages(where="./"),
    # package_dir={"": "./"},
    # cmake_install_dir="./",
    include_package_data=True,
    install_requires=["numpy", "matplotlib", "scipy"],
    cmake_args=["-DBUILD_SHARED_LIBS:BOOL=ON"],
    python_requires=">=3.7",
)
