import numpy as np
from ase import units
from ase.units import Hartree, Bohr


def get_phonon_pol(omega_w, Z_avv, freqs, modes, m_a, cell_cv, eta=0.1e-3):
    # Get phonons at q=0
    Z_vx = Z_avv.swapaxes(0, 1).reshape((3, -1))
    f2_w, D_xw = freqs**2, modes

    alpha_wvv = np.zeros((len(omega_w), 3, 3), dtype=complex)
    m_x = np.repeat(m_a, 3)**0.5
    eta = eta / Hartree
    for f2, D_x in zip(f2_w, D_xw.T):
        if f2 < (1e-3 / Hartree)**2:
            continue
        DM_x = D_x / m_x
        Z_v = np.dot(Z_vx, DM_x)

        alpha_wvv += (np.outer(Z_v, Z_v)[np.newaxis] /
                      ((f2 - omega_w**2) -
                       1j * eta * omega_w)[:, np.newaxis, np.newaxis])

    vol = abs(np.linalg.det(cell_cv))
    alpha_wvv *= 1 / vol

    return alpha_wvv


def phonon_polarizability(bb, Z_avv, m_a, C_NN, cell_cv,
                          symbols, eta=0.1e-3, **overwrite_masses):
    Hartree = units.Hartree
    amuoverme = units._amu / units._me

    for key, value in overwrite_masses.items():
        symbol = key.split('_')[-1]
        for id, symbolid in enumerate(symbols):
            if symbolid == symbol:
                m_a[id] = value

    # Calculate eigenfrequencies
    Minv_NN = np.diag(np.repeat(1 / m_a, 3)**0.5)
    D_NN = np.dot(Minv_NN, np.dot(C_NN, Minv_NN))
    freq2_w, D_xw = np.linalg.eigh(D_NN, UPLO='U')
    s = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
    freq_w = np.sqrt(freq2_w.astype(complex)) * s

    # Make new bb
    bb = dict(bb)
    chiM_qw = bb['chiM_qw']
    chiD_qw = bb['chiD_qw']
    q_abs = bb['q_abs']
    omega_w = bb['omega_w']

    # Get phonons at q=0
    Z_vx = Z_avv.swapaxes(0, 1).reshape((3, -1))
    f2_w = (freq_w / Hartree)**2

    alpha_wvv = np.zeros((len(omega_w), 3, 3), dtype=complex)
    m_x = np.repeat(m_a * amuoverme, 3)**0.5
    gamma = eta / Hartree
    for f2, D_x in zip(f2_w, D_xw.T):
        if f2 < (1e-3 / Hartree)**2:
            continue
        DM_x = D_x / m_x
        Z_v = np.dot(Z_vx, DM_x)

        alpha_wvv += (np.outer(Z_v, Z_v)[np.newaxis] /
                      ((f2 - omega_w**2) -
                       1j * gamma * omega_w)[:, np.newaxis, np.newaxis])

    vol = abs(np.linalg.det(cell_cv)) / units.Bohr**3
    L = np.abs(cell_cv[2, 2] / Bohr)
    alpha_wvv *= 1 / vol * L

    Vm_qw = 2 * np.pi / q_abs[:, None]
    Vd_qw = 2 * np.pi
    newbb = dict(bb)
    chi0Mnew_qw = -(q_abs[:, None])**2 * alpha_wvv[:, 0, 0][np.newaxis]
    chi0Dnew_qw = -(np.ones(len(q_abs))[:, None] *
                    alpha_wvv[:, 2, 2][np.newaxis])
    chi0Mnew_qw += chiM_qw / (1 + Vm_qw * chiM_qw)
    chi0Dnew_qw += chiD_qw / (1 + Vd_qw * chiD_qw)
    chiMnew_qw = chi0Mnew_qw / (1 - Vm_qw * chi0Mnew_qw)
    chiDnew_qw = chi0Dnew_qw / (1 - Vd_qw * chi0Dnew_qw)
    rhoMnew_qz = bb['drhoM_qz']
    rhoDnew_qz = bb['drhoD_qz']

    newbb['chiM_qw'] = chiMnew_qw
    newbb['chiD_qw'] = chiDnew_qw
    newbb['omega_w'] = omega_w
    newbb['q_abs'] = q_abs
    newbb['drhoM_qz'] = rhoMnew_qz
    newbb['drhoD_qz'] = rhoDnew_qz
    newbb['isotropic_q'] = True

    return newbb