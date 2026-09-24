import numpy as np
from assemble import assemble_E
from corner import Jblock
def full_T(m, om, e1, e2=1.0):
    N = len(m.T); k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2); re_ = e1/e2
    E1 = assemble_E(m, k1); E2 = assemble_E(m, k2)
    J = np.array([Jblock(m.n[j], re_, 1.0, 1.0, np.sqrt(re_)) for j in range(N)])
    T = 0.5*(np.eye(8*N) + E2) - 0.5*np.einsum('ixjy,jyz->ixjz', E1.reshape(N, 8, N, 8), J).reshape(8*N, 8*N)
    for j in range(N): T[8*j:8*j+8, 8*j:8*j+8] += 0.5*J[j]
    return T
