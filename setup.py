import os
import subprocess
import setuptools
import unittest
import sys
from pathlib import Path
from enum import Enum, auto

from setuptools.command.install import install
from setuptools.command.test import test

from gramtools.version import package_version


with open("./README.md") as fhandle:
    readme = fhandle.read()

_root_dir = Path(__file__).resolve().parent
cmake_dir = _root_dir / "cmake-build-release"  # cmake target directory

# For pysam to use existing htslib compiled by the backend
os.environ["HTSLIB_LIBRARY_DIR"] = str(cmake_dir / "libgramtools" / "lib")
os.environ["HTSLIB_INCLUDE_DIR"] = str(cmake_dir / "libgramtools" / "include")


class InstallMode(Enum):
    FULL = auto()
    GRAM = auto()


def _build_backend(mode: InstallMode):
    print("Compiling gramtools backend")

    if not os.path.exists(cmake_dir):
        os.mkdir(cmake_dir)
    else:
        print(
            f"Warning: {cmake_dir} already exists; if you want a fresh install, remove it.",
            file=sys.stderr,
        )

    build_type = "REL_WITH_ASSERTS"  # Compiles with -O flag but keeps asserts
    subprocess.call(
        f"CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE={build_type} ..",
        cwd=cmake_dir,
        shell=True,
    )
    target = "gram" if mode is InstallMode.GRAM else "all"

    return_code = subprocess.call(["make", "-j", "4", target], cwd=cmake_dir)
    if return_code != 0:
        print("ERROR: gramtools backend compilation returned: ", return_code)
        exit(-1)


def _test_backend(root_dir):
    print(
        """
    ##########################
    ##Running back end tests##
    ##########################
    """
    )
    test_runner = os.path.join(
        cmake_dir, "libgramtools", "submods", "tests", "test_main"
    )
    test_dir = os.path.join(root_dir, "libgramtools", "submods", "tests")

    return_code = subprocess.call([test_runner], cwd=test_dir)
    if return_code != 0:
        print("ERROR: gramtools backend test runner returned: ", return_code)
        exit(-1)


## Run the front end unit tests: operations making prg, `infer`, integration tests.
# Code adapted from https://stackoverflow.com/q/17001010
def _test_frontend(root_dir):
    print(
        """
    ###########################
    ##Running front end tests##
    ###########################
    """
    )
    # use the default shared TestLoader instance
    test_loader = unittest.defaultTestLoader
    test_suite = test_loader.discover(
        root_dir, pattern="test*.py"
    )  # The pattern used is the default

    # use the basic test runner that outputs to sys.stderr
    # `buffer` suppresses stdout/err if tests pass
    test_runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2, buffer=True)
    tests_results = test_runner.run(test_suite)

    if not tests_results.wasSuccessful():
        print("ERROR: gramtools frontend test runner produced >0 failures")
        exit(-1)


class _InstallCommand(install):
    """
    Command:
        pip3 install -vvv ./gramtools
    """

    def run(self):
        _build_backend(InstallMode.GRAM)
        install.run(self)


class _TestCommand(test):
    """
    Command:
        python3 setup.py test
    """

    def run(self):
        _build_backend(InstallMode.FULL)
        _test_frontend(
            _root_dir
        )  # Front end has integration test relying on backed, so need to build first
        _test_backend(_root_dir)
        print(
            """
        All tests passed.
        """
        )


setuptools.setup(
    name="gramtools",
    version=f"{package_version.__version__}",
    description="Genome inference and variant calling with population reference graphs.",
    url="https://github.com/iqbal-lab-org/gramtools",
    long_description=readme,
    entry_points={"console_scripts": ["gramtools = gramtools.gramtools_main:run"]},
    packages=setuptools.find_packages("."),
    include_package_data=True,
    install_requires=[
        "biopython == 1.76",
        "Cython == 0.29.16",
        "scipy == 1.4.1",
        "pysam == 0.15.4",
        "py-cortex-api == 2.0",
        "cluster_vcf_records >= 0.9.2",
    ],
    test_suite="gramtools.tests",
    cmdclass={"install": _InstallCommand, "test": _TestCommand},
)
