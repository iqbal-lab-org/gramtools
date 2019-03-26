import os
import shutil
import subprocess
import setuptools
import unittest
import sys
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.test import test


with open('./README.md') as fhandle:
    readme = fhandle.read()

_root_dir = os.path.dirname(os.path.realpath(__file__))


def _build_backend(root_dir):
    print('Compiling gramtools backend')
    cmake_dir = os.path.join(root_dir, 'cmake-build-debug')
    try:
        shutil.rmtree(cmake_dir, ignore_errors=False)
    except FileNotFoundError:
        pass
    os.mkdir(cmake_dir)
    subprocess.call('CC=gcc CXX=g++ cmake ..', cwd=cmake_dir, shell=True)

    return_code = subprocess.call(['make'], cwd=cmake_dir)
    if return_code != 0:
        print('ERROR: gramtools backend compilation returned: ', return_code)
        exit(-1)


def _test_backend(root_dir):
    cmake_dir = os.path.join(root_dir, 'cmake-build-debug')
    test_runner = os.path.join(cmake_dir, 'libgramtools', 'tests', 'test_main')
    test_dir = os.path.join(root_dir, 'libgramtools', 'tests')

    return_code = subprocess.call([test_runner], cwd=test_dir)
    if return_code != 0:
        print('ERROR: gramtools backend test runner returned: ', return_code)
        exit(-1)


## Run the front end unit tests: operations making prg, `infer`, `quasimap`.
# Code adapted from https://stackoverflow.com/q/17001010
def _test_frontend(root_dir):
    print("Running gramtools front end unit tests")
    # use the default shared TestLoader instance
    test_loader = unittest.defaultTestLoader
    test_suite = test_loader.discover(root_dir, pattern = "test*.py") # The pattern used is the default

    # use the basic test runner that outputs to sys.stderr
    test_runner = unittest.TextTestRunner(stream = sys.stdout, verbosity = 2)
    tests_results = test_runner.run(test_suite)

    if not tests_results.wasSuccessful():
        print('ERROR: gramtools frontend test runner produced >0 failures')
        exit(-1)


class _InstallCommand(install):
    """
    Command:
        pip3 install -vvv ./gramtools
    """
    def run(self):
        _test_frontend(_root_dir)
        _build_backend(_root_dir)
        _test_backend(_root_dir)
        install.run(self)


class _DevelopCommand(develop):
    """
    Allows for running gramtools as if it were installed (including entry_point), but from source directly.
    Backend is not built on purpose- assuming this will be done using back-end related IDE by developer.
    Tests not ran either- these should be ran manually (eg python3 -m unittest from source root for front-end).

    Command to build develop: (do this inside development dir)
        pip3 install -vvv --editable --process-dependency-links /path/to/gramtools  # wrapper round setuptools develop.

        Note: can also use `python3 /path/to/gramtools/setup.py develop` BUT will not properly run setup.py on packages.

    """
    def run(self):
        develop.run(self)


class _TestCommand(test):
    """
    Command:
        python3 setup.py test
    """
    def run(self):
        _build_backend(_root_dir)
        _test_backend(_root_dir)
        test.run(self)


_package_data = {
    'gramtools': [
        'bin/gram',
        'lib/*',
        'utils/vcf_to_linear_prg.pl'],
}


setuptools.setup(
    name='gramtools',
    version='2.0.0',
    description='Genome inference and variant calling with a reference graph of genetic variation.',
    url='https://github.com/iqbal-lab-org/gramtools',
    long_description=readme,
    entry_points={
        'console_scripts': ['gramtools = gramtools.gramtools:run']
    },
    packages=setuptools.find_packages("."),
    package_data=_package_data,
    include_package_data=True,
    install_requires=[
        'scipy >= 1.0.1',
        'pyvcf >= 0.6.8',
        'py-cortex-api >= 1.0',
        'cluster_vcf_records >= 0.9.2',
    ],
    dependency_links=[
        'https://github.com/iqbal-lab-org/py-cortex-api/tarball/master#egg=py-cortex-api-1.0',
    ],
    test_suite='gramtools.tests',
    cmdclass={
        'install': _InstallCommand,
        'develop': _DevelopCommand,
        'test': _TestCommand,
    })
