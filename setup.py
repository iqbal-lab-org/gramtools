import unittest
import setuptools

# For pysam to use existing htslib compiled by the backend
# os.environ["HTSLIB_LIBRARY_DIR"] = str(cmake_dir / "libgramtools" / "lib")
# os.environ["HTSLIB_INCLUDE_DIR"] = str(cmake_dir / "libgramtools" / "include")


def _test_backend(root_dir):
    print(
        """
    ##########################
    ##Running back end tests##
    ##########################
    """
    )
    test_dir = Path(root_dir) / "libgramtools" / "tests"

    return_code = subprocess.run(
        [str(test_dir / "test_main.bin")], cwd=test_dir
    ).returncode
    if return_code != 0:
        print("ERROR: gramtools backend test runner returned: ", return_code)
        exit(1)


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
        exit(1)


setuptools.setup()
