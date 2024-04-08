from pathlib import Path
import pytest
import os
import qeh
from qeh import interpolate_building_blocks
from qeh import Heterostructure


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

        # WS2
        chi = Path(qeh.__file__).parent / 'chi-data/H-WS2-chi.npz'
        Path(chi.name).symlink_to(chi)

        # WSe2
        chi = Path(qeh.__file__).parent / 'chi-data/H-WSe2-chi.npz'
        Path(chi.name).symlink_to(chi)
        yield tmp_path

    finally:
        os.chdir(cwd)
        print('tmp_path:', tmp_path)


def make_bilayer(matA, matB, thick_A=6.2926, thick_B=6.718, wmax=0):
    # Interpolate to same grid
    interpolate_building_blocks(BBfiles=[matA], BBmotherfile=matB)
    d = (thick_A + thick_B) / 2

    HS = Heterostructure(structure=[matA+'_int', matB+'_int'],
                         d=[d],
                         qmax=None,
                         wmax=wmax,
                         d0=thick_A)
    return HS
