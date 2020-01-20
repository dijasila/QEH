import numpy as np
from ase.units import Bohr, Hartree


# Auxiliary functions
def F(x):
    return x * ((x**2 - 1)**0.5) - np.arccosh(x)


def C(x):
    return x * ((1 - x**2)**0.5) - np.arccos(x)


# Real part of the Polarizability
def P(q_q, w_w, vf, kf):
    prefactor = vf**(-2)
    q_qw = q_q[:, None]
    w_qw = w_w[None, :]
    a = -2 * vf * kf / np.pi
    b = 1 / (4 * np.pi) * \
        (vf * q_qw)**2 / (w_qw**2 - (vf * q_qw)**2)**0.5
    b1 = 1 / (4 * np.pi) * \
        (vf * q_qw)**2 / ((vf * q_qw)**2 - w_qw**2)**0.5
    F1 = F((2 * vf * kf + w_qw) / (vf * q_qw))
    F2 = F((w_qw - 2 * vf * kf) / (vf * q_qw))
    F3 = F((2 * vf * kf - w_qw) / (vf * q_qw))
    C1 = C((2 * vf * kf + w_qw) / (vf * q_qw))
    C2 = C((2 * vf * kf - w_qw) / (vf * q_qw))

    Pol_qw = np.zeros((len(q_q), len(w_w)), complex)
    Pol_qw[:, :] = a

    # Region I
    mask1 = (np.real(w_qw) >= vf * q_qw) * \
            (vf * q_qw + np.real(w_qw) <= 2 * vf * kf)
    Pol_qw += mask1 * b * (F1 - F3)

    # Region II
    mask2 = (np.real(w_qw) >= vf * q_qw) * \
            (vf * q_qw + np.real(w_qw) >= 2 * vf * kf) * \
            (np.real(w_qw) - vf * q_qw <= 2 * vf * kf)
    Pol_qw += mask2 * b * (F1 + 1j * C2)

    # Region III
    mask3 = np.real(w_qw) - vf * q_qw >= 2 * vf * kf
    Pol_qw += mask3 * b * ((F1 - F2) - 1j * np.pi)

    # Region IV
    mask4 = (vf * q_qw >= np.real(w_qw)) * \
            (vf * q_qw + np.real(w_qw) <= 2 * vf * kf)
    Pol_qw += mask4 * 1j * b1 * (F3 - F1)

    # Region V
    mask5 = (vf * q_qw >= np.real(w_qw)) * \
            (vf * q_qw + np.real(w_qw) >= 2 * vf * kf) * \
            (vf * q_qw - 2 * vf * kf <= np.real(w_qw))
    Pol_qw += mask5 * b1 * (C2 - 1j * F1)

    # Region VI
    mask6 = vf * q_qw - 2 * vf * kf >= np.real(w_qw)
    Pol_qw += mask6 * b1 * (C1 + C2)

    Pol_qw *= prefactor

    assert np.allclose(mask1 + mask2 + mask3 + mask4 + mask5 + mask6,
                       1)
    return Pol_qw


# Relaxation time approximation for the Polarizability
def Pgamma(q_q, w_w, tau, vf, kf):
    a = 1j * w_w * tau
    iw_w = w_w + 1j / tau
    P0 = P(q_q, np.array([0j]), vf, kf)
    P1 = P(q_q, iw_w, vf, kf)
    return ((1 - a) * P0 * P1 / (P1 - a * P0))


def GrapheneBB(block, doping, eta):
    Ef = doping / Hartree
    c = 137.0
    vf = c / 300
    kf = Ef / vf
    tau = 1 / (eta / Hartree)
    z = block['z']
    q_q = block['q_abs']
    omega_w = block['omega_w']
    nw = len(omega_w)
    nq = len(q_q)
    chi0M_qw = np.zeros([nq, nw], dtype=complex)
    chiM_qw = np.zeros([nq, nw], dtype=complex)

    V_q = 2 * np.pi / q_q
    chi0M_qw = (Pgamma(q_q, omega_w, tau, vf, kf) +
                (-1.3 / Bohr) * q_q[:, None]**2)
    chiM_qw = chi0M_qw / (1 - chi0M_qw * V_q[:, None])

    # Renormalize monopole density
    drhoM_qz = block['drhoM_qz']

    data = {'isotropic_q': True,
            'q_abs': q_q,
            'omega_w': omega_w,
            'chiM_qw': chiM_qw,
            'chiD_qw': block['chiD_qw'],
            'z': z,
            'drhoM_qz': drhoM_qz,
            'drhoD_qz': block['drhoD_qz']}

    return data
