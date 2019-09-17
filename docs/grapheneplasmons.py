from matplotlib import pyplot as plt
from qeh import make_heterostructure
import numpy as np

layers = ['10graphene+doping=0.1,eta=1e-3']
het = make_heterostructure(layers, frequencies=[1e-5, 0.2, 300],
                           momenta=[0.0001, 0.1, 300])
het.get_plasmon_eigenmodes(filename='bnmodes.npz')

data = np.load('bnmodes.npz')

q_q = data['q_q']
w_w = data['omega_w']
eig_qwl = data['eig']
inveig_qw = - np.sum(1 / eig_qwl, axis=-1).imag

plt.figure(figsize=(3.4, 2.5))
plt.pcolor(q_q, w_w, inveig_qw.T)
ax = plt.gca()
plt.xlabel(r'$q$ (Å$^{-1}$)')
plt.ylabel(r'$\hbar\omega$ (eV)')
plt.tight_layout()
plt.savefig('graphene-multilayer-modes.png', dpi=600)
