from . import version

_, report_dict = version.report()
__version__ = report_dict["version_number"]
