def rotate_bb(block, theta, layer):
    import numpy as np
    from qeh.qeh import interpolate_building_blocks
    from ase.units import Bohr, Hartree
    # Reading old file
    layer = layer
    otherlayers = [layer + '-chi.npz']
    interpolate_building_blocks(BBfiles=otherlayers,
                                q_grid=block['q_abs'] / Bohr,
                                w_grid=block['omega_w'] * Hartree,
                                pad=False)
    block2 = np.load(layer + '_int-chi.npz')
    chiM_qw = block['chiM_qw']
    assert np.allclose(block['q_abs'] - block2['q_abs'], 0.0), (block['q_abs'], block2['q_abs'])
    assert np.allclose(block['z'] - block2['z'], 0.0)
    q_q = block['q_abs']
    omega_w = block['omega_w']

    f1 = np.cos(theta * np.pi / 180)
    f2 = 1 - f1

    chiM_qw = f1 * block['chiM_qw'] + f2 * block2['chiM_qw']
    chiD_qw = f1 * block['chiD_qw'] + f2 * block2['chiD_qw']
    drhoM_qz = f1 * block['drhoM_qz'] + f2 * block2['drhoM_qz']
    drhoD_qz = f1 * block['drhoD_qz'] + f2 * block2['drhoD_qz']

    data = {'isotropic_q': True,
            'q_abs': q_q,
            'omega_w': omega_w,
            'chiM_qw': chiM_qw,
            'chiD_qw': chiD_qw,
            'z': block['z'],
            'drhoM_qz': drhoM_qz,
            'drhoD_qz': drhoD_qz}

    return data
