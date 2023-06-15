import os
from pathlib import Path

import pytest

from qeh import Heterostructure
import qeh


def test_mos2(tmp_path):
    """Calculate static dielectric function for multilayer MoS2
    for 1 to 20 layers
    """
    os.chdir(tmp_path)
    chi = Path(qeh.__file__).parent / 'chi-data/H-MoS2-chi.npz'
    Path(chi.name).symlink_to(chi)

    # positions of maximum:
    q1 = 100000000
    e1 = 0.0

    for n in [1, 2, 3, 4, 5, 10, 20]:
        d = [6.15 for i in range(n - 1)]
        HS = Heterostructure(
            structure=[f'{n}H-MoS2'],  # set up structure
            d=d,  # layer distance array
            include_dipole=True,
            wmax=0,  # only include w=0
            qmax=1,  # q grid up to 1 Ang^{-1}
            d0=6.15)  # width of single layer
        q, w, epsM = HS.get_macroscopic_dielectric_function()

        i = epsM.real.argmax()
        q2 = q[i]
        e2 = epsM.real[i, 0]
        assert q2 < q1  # should decrease as a function of n
        assert e2 > e1  # should increase as a function of n
        q1 = q2
        e1 = e2

    # "Bulk" values:
    assert q1 == pytest.approx(0.118, abs=0.001)
    assert e1 == pytest.approx(12.34, abs=0.01)
