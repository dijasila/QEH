import numpy as np
from conftest import make_bilayer
import pytest


def test_absorption(testdir):
    """Tests absorption spectra for q=0
    """
    HS = make_bilayer('H-MoS2', 'H-WSe2', wmax=10)

    q, w, abs_qw = HS.get_absorption_spectrum()
    previous_absmax = 0  # Needed for flake8...
    previous_wmax = 0
    i = abs_qw[0, :].argmax()
    absmax = abs_qw[0, i]
    wmax = w[i]
    assert absmax == pytest.approx(313.324, abs=0.001)
    assert wmax == pytest.approx(2.958, abs=0.001)

    # Check so that all values are positive
    # but could be small negative parts due to numerics
    assert np.all(abs_qw >= -1e-6)


@pytest.mark.xfail
def test_absorption_q(testdir):
    """Tests absorption spectra finite q.
    This test fails, some peaks at finit q have
    higher intensity and peak position does not shift to
    higher frequency. Is this correct?
    """
    HS = make_bilayer('H-MoS2', 'H-WSe2', wmax=10)

    q, w, abs_qw = HS.get_absorption_spectrum()
    previous_absmax = 0  # Needed for flake8...
    previous_wmax = 0
    for iq in range(10):
        i = abs_qw[iq, :].argmax()
        absmax = abs_qw[iq, i]
        wmax = w[i]
        if iq == 0:
            # Test values gamma point
            assert absmax == pytest.approx(313.324, abs=0.001)
            assert wmax == pytest.approx(2.958, abs=0.001)
        else:
            # Test so that perak position shifts to higher energy
            # and max value decreases
            assert absmax < previous_absmax
            assert wmax > previous_wmax
        previous_absmax = absmax
        previous_wmax = wmax

    # Check so that all values are positive
    # but could be small negative parts due to numerics
    assert np.all(abs_qw >= -1e-6)
