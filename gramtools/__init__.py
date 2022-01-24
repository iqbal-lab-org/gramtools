from pathlib import Path

from . import version

_, report_dict = version.report()
__version__ = report_dict["version_number"]

py_cortex_api_message = """
To run discovery, you need the following:
    - 'pip install py-cortex-api==2.2.0'
    - R and perl on your $PATH
        """


# Executable/library locations
_base_install_path = Path(__file__).resolve().parent
gramtools_exec_fpath = str(_base_install_path / "bin" / "gram")
gramtools_lib_fpath = str(_base_install_path / "lib")
