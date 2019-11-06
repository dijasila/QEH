import numpy as np
from numpy.lib.scimath import sqrt as csqrt
from ase.units import Hartree


# Chi0 for zero temperature
def chi0T0(q_qwm, w_qwm, me, mup_qwm):
    kf = np.sqrt(2 * me * mup_qwm)
    N = kf**2 / (2 * np.pi)
    vf = kf / me
    z = q_qwm / (2 * kf)
    u = w_qwm / (q_qwm * vf)
    G = N / (me * z * (vf**2))

    # Condition 1
    mask1 = (np.absolute(z - u.real) >= 1)
    chi0T0_qwm = (mask1 * (-(z - u.real) / np.absolute(z - u.real)) *
                  csqrt((z - u)**2 - 1))

    # Condition 2
    mask2 = (np.absolute(z + u.real) >= 1)
    chi0T0_qwm -= (mask2 * ((z + u.real) / np.absolute(z + u.real)) *
                   csqrt((z + u)**2 - 1))

    # Condition 3
    mask3 = (np.absolute(z - u.real) < 1)
    chi0T0_qwm += 1j * mask3 * csqrt(1 - (z - u)**2)

    # Condition 4
    mask4 = (np.absolute(z + u.real) < 1)
    chi0T0_qwm -= 1j * mask4 * csqrt(1 - (z + u)**2)

    chi0T0_qwm += 2 * z
    chi0T0_qwm *= G
    return -chi0T0_qwm


# Sum Argument
def arg(q, w, me, T, mu, mup):
    argument = (chi0T0(q, w, me, mup) /
                (4 * T * (np.cosh((mu - mup) / (2 * T)))**2))
    return argument


# Polarizability
def P(q_q, w_w, me, T, mu, mupmax, N):
    mup_m = np.linspace(10**(-5), mupmax, N)
    q_qwm = q_q[:, None, None]
    w_qwm = w_w[None, :, None]
    mup_qwm = mup_m[None, None, :]
    return np.trapz(arg(q_qwm, w_qwm, me, T, mu, mup_qwm), mup_qwm, axis=2)


# Polarizability in the relaxation time approximation
def Pgamma(q_q, w_w, me=None, efermi=None, T=0.0,
           mupmax=None, N=1000, gamma=1e-4):
    assert efermi is not None, print('You have to set a fermi energy!')
    assert me is not None, \
        print('You have to set an effective electron mass!')

    gamma = gamma / Hartree
    efermi = efermi / Hartree
    T = T / Hartree
    a = 1j * gamma / w_w
    iw_w = w_w + 1j * gamma

    if T / efermi > 1e-7:
        # Temperature dependent chemical potential
        mu = T * np.log(np.exp(efermi / T) - 1)
        if mupmax is None:
            mupmax = 20 * T + mu * (mu > 0)
        P0 = P(q_q, np.array([0j]), me, T, mu, mupmax, N)
        P1 = P(q_q, iw_w, me, T, mu, mupmax, N)
    else:
        mu = efermi
        T = 0
        P0 = chi0T0(q_q[:, None], np.array([0j])[None, :], me, mu)
        P1 = chi0T0(q_q[:, None], iw_w[None, :], me, mu)

    return ((1 + a) * P1 / (1 + a * P1 / P0))


def dopedsemiconductor(block, mex, mey, doping,
                       temperature, eta, theta=0):
    # Reading old file
    chiM_qw = block['chiM_qw']
    qgrid_q = block['q_abs']
    omega_w = block['omega_w']
    me = np.sqrt(mex * mey)
    trans = [(mey / mex)**(1 / 4), (mex / mey)**(1 / 4)]
    qx = np.cos(theta / 180 * np.pi) * qgrid_q
    qy = np.sin(theta / 180 * np.pi) * qgrid_q
    q_q = np.sqrt((trans[0] * qx)**2 + (trans[1] * qy)**2)
    V_q = 2 * np.pi / qgrid_q
    chi0Mnew_qw = Pgamma(q_q, omega_w, me=me,
                         efermi=doping, T=temperature,
                         gamma=eta)
    chi0Mnew_qw += chiM_qw / (1 + V_q[:, None] * chiM_qw)
    dopedchiM_qw = chi0Mnew_qw / (1 - chi0Mnew_qw * V_q[:, None])

    data = {'isotropic_q': True,
            'q_abs': qgrid_q,
            'omega_w': omega_w,
            'chiM_qw': dopedchiM_qw,
            'chiD_qw': block['chiD_qw'],
            'z': block['z'],
            'drhoM_qz': block['drhoM_qz'],
            'drhoD_qz': block['drhoD_qz']}

    return data
