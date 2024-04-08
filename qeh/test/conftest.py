from pathlib import Path
import pytest
import os
import qeh


@pytest.fixture
def testdir(tmp_path):
    """Create and go to temporary directory."""

    cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        # Add softlinks to building blocks used in tests
        # MoS2
        chi = Path(qeh.__file__).parent / 'chi-data/H-MoS2-chi.npz'
        Path(chi.name).symlink_to(chi)

        #WS2
        chi = Path(qeh.__file__).parent / 'chi-data/H-WS2-chi.npz'
        Path(chi.name).symlink_to(chi)
        yield tmp_path
    finally:
        os.chdir(cwd)
        print('tmp_path:', tmp_path)
