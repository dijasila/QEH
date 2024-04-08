from pathlib import Path
import pytest
import os


@pytest.fixture
def testdir(tmp_path):
    """Create and go to temporary directory."""
    cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        yield tmp_path
    finally:
        os.chdir(cwd)
        print('tmp_path:', tmp_path)
