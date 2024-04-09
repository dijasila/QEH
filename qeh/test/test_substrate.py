import numpy as np
from qeh import Heterostructure
import pytest


def mock_substrate(d):
    ei = 2.4
    f1 = 0.7514
    f2 = 0.1503
    f3 = 0.6011
    w1 = 0.055
    w2 = 0.098
    w3 = 0.140
    w = np.linspace(start=0, stop=1.3, num=2000)
    A = f1 * w1**2 / (w1**2+w**2)
    B = f2 * w2**2 / (w2**2+w**2)
    C = f3 * w3**2 / (w3**2+w**2)
    eps_w = ei + A + B + C
    substrate = {'eps': eps_w, 'omega': w, 'd': [d], 'isotropic': True}
    return substrate


def test_substrate(testdir):

    # First test isolated layer
    d_MoS2 = 6.15
    hl_array = np.array([1, 0])
    el_array = np.array([1, 0])
    HS = Heterostructure(structure=['1H-MoS2'],
                         d=[],
                         wmax=0,
                         d0=d_MoS2)
    ee, ev = HS.get_exciton_binding_energies(eff_mass=0.27,
                                             e_distr=el_array,
                                             h_distr=hl_array)
    Eb_isolated_layer = -ee[0]
    assert Eb_isolated_layer == pytest.approx(0.592, abs=0.001)

    # Now add substrate
    substrate = mock_substrate(d_MoS2 / 2)
    HS_sub = Heterostructure(structure=['1H-MoS2'],
                             d=[],
                             wmax=0,
                             d0=d_MoS2,
                             substrate=substrate)
    ee, ev = HS_sub.get_exciton_binding_energies(eff_mass=0.27,
                                                 e_distr=el_array,
                                                 h_distr=hl_array)
    Eb_sub = -ee[0]
    assert Eb_sub == pytest.approx(0.354, abs=0.001)

    # Increase substrate distance and verify that exciton BE increases
    substrate = mock_substrate(d_MoS2)
    HS_sub = Heterostructure(structure=['1H-MoS2'],
                             d=[],
                             wmax=0,
                             d0=d_MoS2,
                             substrate=substrate)
    ee, ev = HS_sub.get_exciton_binding_energies(eff_mass=0.27,
                                                 e_distr=el_array,
                                                 h_distr=hl_array)
    Eb_sub2 = -ee[0]
    assert Eb_sub2 == pytest.approx(0.4704, abs=0.001)
