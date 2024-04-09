import numpy as np
from conftest import make_bilayer
import pytest


def test_eels(testdir):
    """Calculate eels spectra
    """
    HS = make_bilayer('H-MoS2', 'H-WSe2', wmax=10)

    q, w, eels_qw = HS.get_eels()
    previous_eelsmax = 0  # Needed for flake8...
    previous_wmax = 0

    for iq in range(10):
        i = eels_qw[iq, :].argmax()
        eelsmax = eels_qw[iq, i]
        wmax = w[i]
        if iq == 0:
            # Test values gamma point
            assert eelsmax == pytest.approx(156.613, abs=0.001)
            assert wmax == pytest.approx(2.958, abs=0.001)
        else:
            # Test so that perak position shifts to higher energy
            # and max value decreases
            assert eelsmax < previous_eelsmax
            assert wmax > previous_wmax
        previous_eelsmax = eelsmax
        previous_wmax = wmax

    # Check so that all values are positive
    # but could be small negative parts due to numerics
    assert np.all(eels_qw >= -1e-6)
