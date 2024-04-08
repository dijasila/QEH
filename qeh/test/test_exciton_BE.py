import numpy as np
from conftest import make_bilayer


def test_exciton_BE(testdir):
    """Calculate W(q)
    """
    HS = make_bilayer('H-MoS2', 'H-WSe2')

    hl_array = np.array([1., 0., 0, 0.])
    el_array = np.array([0., 0., 1., 0.])
    inter_mass = 0.244
    ee, ev = HS.get_exciton_binding_energies(eff_mass=inter_mass,
                                             e_distr=el_array,
                                             h_distr=hl_array)
    
    assert np.allclose(-ee[0], 0.28068, atol=1e-4)
