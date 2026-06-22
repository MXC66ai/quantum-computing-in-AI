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
