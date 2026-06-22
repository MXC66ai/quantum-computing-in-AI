import numpy as np
from qiskit_dynamics import Solver, LindbladModel
from qiskit.quantum_info import Operator

# adjacency matrix for 5-node graph (HIF-1α + 3 genes + PolII)
M = np.array([[ -0.2,  0.1, 0.05, 0.05, 0.0],
              [ 0.0, -0.1, 0.0,  0.0,  0.1],
              [ 0.0, 0.0, -0.08, 0.0,  0.08],
              [ 0.0, 0.0, 0.0, -0.07, 0.07],
              [ 0.0, 0.0, 0.0, 0.0,  0.0]])

# gravitational scaling factor
gamma_g = 1.0           # 1g
# gamma_g = 0.001        # μg

M_mu = gamma_g * M       # apply scaling for μg condition

# Lindblad operators (decoherence)
L = []                  # list of jump operators (optional, here set empty)

ham = -1j*(M_mu - M_mu.T)    # coherent part
model = LindbladModel(hamiltonian=ham, dissipators=L)
solver = Solver(model)

t_grid = np.linspace(0, 600, 61)   # 0–600 s, 10 s step
psi0 = np.zeros(5); psi0[0] = 1.0   # start with HIF‑1α only
result = solver.solve(t_grid, psi0)

# Extract probabilities of each gene ON (state = gene)
P_Tef = result.y[:,1]   # index 1 = Tef
P_Sst = result.y[:,2]
P_Oas = result.y[:,3]

print("Final ON probabilities (μg):", P_Tef[-1], P_Sst[-1], P_Oas[-1])
