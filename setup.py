import os
import shutil
import subprocess
import setuptools
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.test import test


with open('./README.md') as fhandle:
    readme = fhandle.read()


class LibGramtools:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.cmake_dir = os.path.join(self.root_dir, 'cmake-build-debug')

    def build(self):
        print('Compiling gramtools backend')
        try:
            shutil.rmtree(self.cmake_dir, ignore_errors=False)
        except FileNotFoundError:
            pass
        subprocess.call(['mkdir', self.cmake_dir])

        subprocess.call('CC=gcc CXX=g++ cmake ..', cwd=self.cmake_dir, shell=True)
        subprocess.call(['make'], cwd=self.cmake_dir)

    def test(self):
        test_runner = os.path.join(self.cmake_dir,
                                   'libgramtools',
                                   'tests',
                                   'test_main')
        test_dir = os.path.join(self.root_dir,
                                'libgramtools',
                                'tests')

        return_code = subprocess.call([test_runner], cwd=test_dir)
        if return_code != 0:
            print('ERROR: libgramtools test runner returned: ', return_code)
            exit(-1)


class InstallCommand(install):
    """pip3 install -vvv ./gramtools"""
    def run(self):
        libgramtools = LibGramtools()
        libgramtools.build()
        libgramtools.test()
        install.run(self)


class DevelopCommand(develop):
    """pip3 install -vvv --editable ./gramtools"""
    def run(self):
        libgramtools = LibGramtools()
        libgramtools.build()
        libgramtools.test()
        develop.run(self)


class TestCommand(test):
    """python3 setup.py test"""
    def run(self):
        libgramtools = LibGramtools()
        libgramtools.build()
        libgramtools.test()
        test.run(self)


package_data = {
    'gramtools': [
        'bin/gram',
        'lib/*',
        'utils/vcf_to_linear_prg.pl'],
}


setuptools.setup(
    name='gramtools',
    version='2.0',
    description='Genome inference and variant calling with a reference graph of genetic variation.',
    url='https://github.com/iqbal-lab-org/gramtools',
    long_description=readme,
    entry_points={
        'console_scripts': ['gramtools = gramtools.gramtools:run']
    },
    packages=setuptools.find_packages("."),
    package_data=package_data,
    include_package_data=True,
    test_suite='gramtools.tests',
    cmdclass={
        'install': InstallCommand,
        'develop': DevelopCommand,
        'test': TestCommand,
    })
