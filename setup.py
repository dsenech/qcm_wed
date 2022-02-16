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
