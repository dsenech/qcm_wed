import sys
import os
from setuptools import find_packages
from skbuild import setup

#ne pas importer numpy, et pas de try/except sur skbuild, mais plutot mentionner son installation dans le readme

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
    license="GPL",
    packages=find_packages(),
    # package_dir={"": "."},
    cmake_install_dir="./",
    include_package_data=True,
    install_requires=["numpy", "matplotlib", "scipy"],
    python_requires=">=3.7",
)
