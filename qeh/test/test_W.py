from qeh import interpolate_building_blocks
from qeh import Heterostructure
import numpy as np


def test_W(testdir):
    """Calculate W(q)
    """

    # Interpolate to same grid
    interpolate_building_blocks(BBfiles=['H-WS2'], BBmotherfile='H-MoS2')

    thick_MoS2 = 6.2926
    thick_WSe2 = 6.718

    d_MoS2_WSe2 = (thick_MoS2 + thick_WSe2) / 2

    HS = Heterostructure(structure=['H-MoS2_int', 'H-WS2_int'],
                         d=[d_MoS2_WSe2],
                         qmax=None,
                         wmax=0,
                         d0=thick_WSe2)

    hl_array = np.array([1., 0., 0, 0.])
    el_array = np.array([1., 0., 0., 0.])

    # intralayer exciton potential
    q, Wintra_q, _ = HS.get_exciton_screened_potential(el_array, hl_array)

    # average intralayer screened potential
    Wintra_qw = HS.get_screened_potential(layer=0)

    # Exciton potential should be minus average potential
    assert np.allclose(Wintra_qw[:, 0], -Wintra_q)

    # Check some numbers as well...
    assert np.allclose([Wintra_q[0], Wintra_q[50], Wintra_q[86]],
                       [-626451.0099793368, -0.8781853200860702,
                        -0.3942195574531765], atol=1e-5)

    # Getting the interlayer exciton screened potential
    # on real grid
    hl_array = np.array([0., 0., 1., 0.])
    el_array = np.array([1., 0., 0., 0.])

    r, W_r = HS.get_exciton_screened_potential_r(
        r_array=np.linspace(1e-1, 30, 1000),
        e_distr=el_array,
        h_distr=hl_array)
    assert np.allclose([W_r[0], W_r[-1]], [-0.0190095, -0.0120897],
                       atol=1e-5)
