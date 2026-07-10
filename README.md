<img width="611" height="135" alt="image" src="https://github.com/user-attachments/assets/bdfaa70f-b336-4356-9ada-caca383cd89f" />

example 1:
from qiskit_nature.drivers import PySCFDriver
from qiskit_nature.transformers import FreezeCoreTransformer
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import ParityMapper

driver = PySCFDriver(atom=geom, unit='Angstrom', charge=0, spin=0)
problem = driver.run()
transformer = FreezeCoreTransformer()
electronic_problem = transformer.transform(problem)
converter = QubitConverter(ParityMapper(), two_qubit_reduction=True)
qubit_op = converter.convert(electronic_problem.second_q_ops()[0])

example 2_Sample VQE script:
from qiskit import Aer, IBMQ
from qiskit.utils import QuantumInstance
from qiskit.algorithms import VQE
from qiskit.circuit.library import EfficientSU2
from qiskit.algorithms.optimizers import SLSQP

backend = Aer.get_backend('qasm_simulator')
quantum_instance = QuantumInstance(backend, shots=5000)

ansatz = EfficientSU2(num_qubits=qubit_op.num_qubits, entanglement='full')
optimizer = SLSQP(maxiter=200)

vqe = VQE(ansatz, optimizer, quantum_instance=quantum_instance)
energy = vqe.compute_minimum_eigenvalue(qubit_op).eigenvalue.real
print(f"VQE ground‑state energy = {energy:.6f} Ha")


Example QSW code :
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



Sample Pandoc Conversion Commands:
# Convert markdown → PDF (LaTeX engine = xelatex)
pandoc Quantum_Simulation_HIF1a_Protocol.md \
    -o Quantum_Simulation_HIF1a_Protocol.pdf \
    --pdf-engine=xelatex \
    -V geometry:margin=1in \
    -V fontsize=12pt

# Convert markdown → DOCX
pandoc Quantum_Simulation_HIF1a_Protocol.md \
    -o Quantum_Simulation_HIF1a_Protocol.docx


P.S.需要安装必要软件（`pandoc`, `texlive-xetex`, `python`, `qiskit`）
