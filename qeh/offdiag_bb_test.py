import numpy as np
from gpaw.response.df import DielectricFunction
from gpaw.response.qeh import BuildingBlock
from gpaw import GPAW, PW, FermiDirac
from qeh import Heterostructure, interpolate_building_blocks
from matplotlib import pyplot as plt
import os
from pathlib import Path
import qeh
from ase.units import Bohr


def dielectric(calc, domega, omega2, rate=0.0):
    diel = DielectricFunction(calc=calc,
                              frequencies={'type': 'nonlinear',
                                           'omegamax': 10,
                                           'domega0': domega,
                                           'omega2': omega2},
                              nblocks=1,
                              ecut=10,
                              rate=rate,
                              truncation='2D',
                              txt='df.txt')
    return diel


def IBiTe_gs():
    # janus material. material parameters obtained from c2db.
    from ase.atoms import Atoms
    IBiTe_positions = np.array([[0, 2.552, 7.802],
                                [0, 0, 9.872],
                                [2.210, 1.276, 11.575]])
    IBiTe = Atoms('IBiTe', positions=IBiTe_positions)
    IBiTe.pbc = [True, True, False]
    cell = np.array([[4.4219, 0, 0.0, ],
                     [-2.211, 3.829, 0.0],
                     [0.0, 0.0, 19.5]])
    IBiTe.cell = cell
    calc = GPAW(mode=PW(200),
                xc='LDA',
                occupations=FermiDirac(0.01),
                kpts={'size': (6, 6, 1), 'gamma': True},
                txt=None)
    IBiTe.calc = calc
    IBiTe.get_potential_energy()
    calc.write('IBITe.gpw', mode='all')
    q_cs = np.array([[0, 0, 0], [1/6, 0, 0], [2/6, 0, 0], [3/6, 0, 0]])
    rcell_cv = 2 * np.pi * np.linalg.inv(calc.wfs.gd.cell_cv).T
    q_vs = np.dot(q_cs, rcell_cv)
    q_abs=[np.linalg.norm(q_v) / Bohr for q_v in q_vs]
    return q_cs, q_abs


def make_HS(structure, off_diag):
    # XXX Monolayer results are very sensitive to d0
    # GPAW results are retained if d0 is chosen as
    # lattice constant in z-direction.
    # investigate the reason for this...
    d = [7]
    # Including off-diagonal elements (True by default)
    HS = Heterostructure(
        structure=structure,  # set up structure
        d=d,  # layer distance array
        include_dipole=True,
        include_off_diagonal=off_diag,
        wmax=0,  # only include w=0
        qmax=1,  # q grid up to 1 Ang^{-1}
        d0=19.5)  # width of single layer
    return HS


def test_off_diagonal_chi(tmp_path):
    # Calculating gs and buildingblock
    # XXX Move to fixture as we create more tests
    os.chdir(tmp_path)
    chi = Path(qeh.__file__).parent / 'chi-data/H-MoS2-chi.npz'
    Path(chi.name).symlink_to(chi)

    q_cs, q_abs = IBiTe_gs()
    df = dielectric('IBiTe.gpw', 0.1, 0.5)
    bb = BuildingBlock('IBiTe', df)
    bb.calculate_building_block()

    # Compare monolayer with gpaw
    monolayer = make_HS(['IBiTe'], False)
    q, w, epsM_mono = monolayer.get_macroscopic_dielectric_function()

    epsM_gpaw = []
    tested_qs = 0
    for iq, q_c in enumerate(q_cs):
        _, epsM_q = df.get_dielectric_function(q_c=q_c)
        epsM_gpaw.append(epsM_q[0])
        for iq2 in range(len(q)):
            if np.isclose(q[iq2], q_abs[iq], atol=5e-5):
                assert np.isclose(epsM_q[0], epsM_mono[iq2], atol=0.03)
                tested_qs += 1
    assert tested_qs == 3
    if False:
        # XXX possibility to plot, remove when satisfied...
        epsM_gpaw = np.array(epsM_gpaw)
        plt.plot(q, epsM_mono.real, '-*', label='qeh real')
        plt.plot(q, epsM_mono.imag, '-*', label='qeh imag')
        plt.plot(q_abs, epsM_gpaw.real, '-+', label='gpaw real')
        plt.plot(q_abs, epsM_gpaw.imag, '-+', label='gpaw imag')
        plt.legend()
        plt.title('Macroscopic dielectric function')
        plt.show()

    # Including off-diag bb
    HS = make_HS(['2IBiTe'], True)
    q, w, epsM = HS.get_macroscopic_dielectric_function()

    # excluding off-diag bb
    HS2 = make_HS(['2IBiTe'], False)
    q2, w2, epsM2 = HS2.get_macroscopic_dielectric_function()

    # off diag building block should make a small but
    # finite difference
    assert not np.allclose(epsM, epsM2, atol=1e-1)
    assert np.allclose(epsM, epsM2, atol=5e-1)

    # Test actual values
    # XXX It does not seem correct to have imaginary parts here?
    # However, this is the case also when off_diag = False, so
    # maybe it is the building block which is wrong?
    expected_epsM = [[1.00134782 - 2.86604242e-06j],
                     [1.92047859 - 4.27676734e-03j],
                     [2.79431902 - 1.17221028e-02j],
                     [3.71008573 - 2.00242298e-02j],
                     [4.75955272 - 2.49296881e-02j],
                     [6.06830419 - 1.72620056e-02j],
                     [7.84803204 + 2.52031425e-02j],
                     [10.52641274 + 1.65458907e-01j],
                     [15.16435212 + 6.30124718e-01j],
                     [25.34195424 + 2.65514592e+00j],
                     [5.401448 + 4.28917276e-02j],
                     [2.27880957 + 4.93372479e-02j]]
    assert np.allclose(epsM, expected_epsM, atol=1e-4)

    # Test combination Janus + non Janus material
    # dielectric function should depend on stacking
    interpolate_building_blocks(BBfiles=['IBiTe', 'H-MoS2'],
                                BBmotherfile='H-MoS2')

    HS3 = make_HS(['IBITe_int', 'H-MoS2_int'], False)
    q3, w3, epsM3 = HS3.get_macroscopic_dielectric_function()
    HS4 = make_HS(['H-MoS2_int', 'IBITe_int'], False)
    q4, w4, epsM4 = HS4.get_macroscopic_dielectric_function()

    # Check that maximum values of epsM are correct, and different for the two
    # stackings. Note this is mainly due to density and not so much due to off-
    # diag bb
    assert np.isclose(np.amax(epsM3.real), 12.2796, atol=1e-3)
    assert np.isclose(np.amax(epsM4.real), 11.7102, atol=1e-3)
