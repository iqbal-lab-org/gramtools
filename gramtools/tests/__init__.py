from pathlib import Path

from gramtools import gramtools_main

base_dir = Path(gramtools_main.__file__).resolve().parent
data_dir = base_dir / "tests" / "integration_test_data"

gramtools_main._setup_parser()
