import time
import logging
import collections
import os
import json

from gramtools import version
from .paths import ProjectPaths

log = logging.getLogger("gramtools")


def with_report(f):
    """
    Decorator to add logging and reporting to gramtools procedures
    To signal that something went wrong in decorated function call, that function needs to raise an Exception.
    """

    def reportify(report, action: str, command_paths: ProjectPaths, *args):
        success, error, error_message = True, None, None
        timer_start = time.time()

        try:
            original_result = f(report, action, command_paths, *args)
        except Exception as e:
            success = False
            error = e
        timer_end = time.time()

        report["success"] = success
        process_report = collections.OrderedDict(
            {
                "success": success,
                "error_message": str(error),
                "run_time": int(timer_end) - int(timer_start),
            }
        )

        # The condition below accounts for the called function reporting as well.
        if action not in report:
            report["processes"][action] = process_report
        else:
            report["processes"][action] = {
                **process_report,
                **report[action],
            }  # Place success status at very top

        if not success:
            log.error(f"{error.__class__.__name__} raised with message: {error}")
            log.error(
                f"Unsuccessful {action}. " f"Process reported to {command_paths.report}"
            )
            _save_report(report, command_paths)
            exit(1)

        log.debug(f"Ran {action} in: {timer_end - timer_start} seconds")
        return original_result

    return reportify


def _save_report(report, command_paths: ProjectPaths, command_hash_paths=None):
    log.debug("Saving command report:\n%s", command_paths.report)

    end_time = str(time.time()).split(".")[0]
    _, report_dict = version.report()
    current_working_directory = os.getcwd()

    start_time = report.pop("start_time")
    report.update(
        collections.OrderedDict(
            [
                ("total_runtime", int(end_time) - int(start_time)),
                ("current_working_directory", current_working_directory),
                ("paths", command_paths.dict()),
                ("path_hashes", command_hash_paths),
                ("version_report", report_dict),
            ]
        )
    )

    with open(command_paths.report, "w") as fhandle:
        json.dump(report, fhandle, indent=4)


def new_report():
    start_time = str(time.time()).split(".")[0]
    rep = collections.OrderedDict(
        {
            "success": "",
            "processes": collections.OrderedDict(),
            "start_time": start_time,
        }
    )
    return rep
