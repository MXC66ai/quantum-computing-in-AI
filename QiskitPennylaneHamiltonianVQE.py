#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ====================  IMPORTS ====================
import numpy as np, pandas as pd
from qiskit import QuantumCircuit, Aer, transpile
from qiskit.opflow import I, X, Y, Z, PauliSumOp
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import SPSA
from qiskit.primitives import Sampler
from qiskit.circuit.library import EfficientSU2
from qiskit.utils import algorithm_globals
from qiskit.providers.fake_provider import FakeMelbourne
import matplotlib.pyplot as plt

# ====================  1. DATA LOADING ====================
node_list = ['CXCL12','CXCL8','CXCR4','CXCR2','ACKR3']
node_idx = {n:i for i,n in enumerate(node_list)}

# ---- 1.1 边信息（示例） ----
edge_df = pd.read_csv('data/edge_weights.csv')    # 见前文格式
edge_df['J_raw'] = edge_df['weight'] / edge_df['Kd']
J_max = edge_df['J_raw'].max()
edge_df['J_norm'] = edge_df['J_raw'] / J_max
edge_list = [(node_idx[r.ligand], node_idx[r.receptor], r.J_norm)
             for _,r in edge_df.iterrows()]

# ---- 1.2 节点信息 ----
node_exp = pd.read_csv('data/node_expr.csv')      # columns: node, expr
expr_norm = (node_exp['expr']-node_exp['expr'].min()) / \
            (node_exp['expr'].max()-node_exp['expr'].min())
h_angle = (expr_norm - 0.5) * np.pi                # [-π/2, +π/2]
h_dict = {node_idx[row.node]: angle for row,angle in zip(node_exp.itertuples(),h_angle)}

# 横向场 (可根据 IL‑6 ELISA 决定)
gamma_dict = {i: 0.2 for i in range(len(node_list))}

# ==================== 2. HAMILTONIAN ====================
def make_tfim_hamiltonian(num_qubits, h_dict, J_edges, gamma_dict):
    # 1) Z term
    H = sum(h_dict[i] * (Z if i==0 else I^(i-1) ^ Z ^ I^(num_qubits-i-1))
            for i in range(num_qubits))
    # 2) X term
    H += sum(gamma_dict[i] * (X if i==0 else I^(i-1) ^ X ^ I^(num_qubits-i-1))
             for i in range(num_qubits))
    # 3) XX+YY coupling
    for (i,j,J) in J_edges:
        term_xx = J * (X if i==0 else I^(i-1) ^ X ^ I^(num_qubits-i-1))
        term_xx = term_xx.compose( X if j==0 else I^(j-1) ^ X ^ I^(num_qubits-j-1) )
        term_yy = J * (Y if i==0 else I^(i-1) ^ Y ^ I^(num_qubits-i-1))
        term_yy = term_yy.compose( Y if j==0 else I^(j-1) ^ Y ^ I^(num_qubits-j-1) )
        H += term_xx + term_yy
    return PauliSumOp.from_operator(H.to_matrix())

num_qubits = len(node_list)
H_op = make_tfim_hamiltonian(num_qubits, h_dict, edge_list, gamma_dict)

# ==================== 3. VQE ====================
algorithm_globals.random_seed = 42
ansatz = EfficientSU2(num_qubits, entanglement='full', reps=2)
optimizer = SPSA(maxiter=150)

# 使用 Qiskit Aer simulator（可以换成真实后端）
backend = Aer.get_backend('aer_simulator_statevector')
sampler = Sampler(backend=backend)

vqe = VQE(ansatz=ansatz,
          optimizer=optimizer,
          sampler=sampler)

result = vqe.compute_minimum_eigenvalue(operator=H_op)
E0 = result.eigenvalue.real
print(f"[VQE] 基态能量 = {E0:.6f}")

# ------------------- 计算磁化 -------------------
opt_params = result.optimal_point
circuit_opt = ansatz.assign_parameters(opt_params)
circuit_opt.measure_all()
job = backend.run(transpile(circuit_opt, backend))
counts = job.result().get_counts()
# 直接用状态向量更简洁（这里用 statevector）
state = backend.run(circuit_opt).result().get_statevector()
mag = {}
for i,name in enumerate(node_list):
    Z_i = (Z if i==0 else I^(i-1) ^ Z ^ I^(num_qubits-i-1)).to_matrix()
    mag[name] = np.real(state.expectation_value(Z_i))
print("节点磁化 (⟨σ_z⟩):", mag)

M_cxcl = np.mean([mag[n] for n in ['CXCL12','CXCL8']])
print(f"CXCL 总磁化 M = {M_cxcl:.4f}")

# ==================== 4. TIME‑EVOLUTION (TROTTER) ====================
def trotter_step(circ, h_dict, J_edges, gamma_dict, dt):
    for i, h in h_dict.items():
        circ.rz(2*h*dt, i)
    for i, g in gamma_dict.items():
        circ.rx(2*g*dt, i)
    for (i,j,J) in J_edges:
        circ.rxx(2*J*dt, i, j)
        circ.ryy(2*J*dt, i, j)
    return circ

def build_evolution(num_qubits, h_dict, J_edges, gamma_dict,
                   total_t=10.0, n_steps=20):
    dt = total_t / n_steps
    circ = QuantumCircuit(num_qubits)
    circ.h(range(num_qubits))       # 产生均匀叠加
    for _ in range(n_steps):
        circ = trotter_step(circ, h_dict, J_edges, gamma_dict, dt)
    circ.measure_all()
    return circ

evo_circ = build_evolution(num_qubits, h_dict, edge_list,
                            gamma_dict, total_t=30.0, n_steps=30)
print(evo_circ.draw(output='mpl'))

# 在真实硬件上跑（示例）：
# backend_real = provider.get_backend('ibmq_16_melbourne')
# sampler_real = Sampler(backend=backend_real, options={"shots":8192, "resilience_level":1})
# result_evo = sampler_real.run(evo_circ).result()
# ...

# ==================== 5. QUANTUM WALK (LAPLACIAN) ====================
def laplacian_hamiltonian(num_qubits, W):
    deg = np.sum(W, axis=1)
    H = sum(deg[i] * (Z if i==0 else I^(i-1) ^ Z ^ I^(num_qubits-i-1))
            for i in range(num_qubits))
    for i in range(num_qubits):
        for j in range(i+1, num_qubits):
            J = W[i,j]
            if J==0: continue
            term_xx = -J * (X if i==0 else I^(i-1) ^ X ^ I^(num_qubits-i-1))
            term_xx = term_xx.compose( X if j==0 else I^(j-1) ^ X ^ I^(num_qubits-j-1))
            term_yy = -J * (Y if i==0 else I^(i-1) ^ Y ^ I^(num_qubits-i-1))
            term_yy = term_yy.compose( Y if j==0 else I^(j-1) ^ Y ^ I^(num_qubits-j-1))
            H += term_xx + term_yy
    return PauliSumOp.from_operator(H.to_matrix())

# 生成对称化 W（只取上三角）
W_mat = np.zeros((num_qubits, num_qubits))
for i,j,J in edge_list:
    W_mat[i,j] = J
    W_mat[j,i] = J

L_op = laplacian_hamiltonian(num_qubits, W_mat)

# 用 VQE 求最小非零特征值 λ1
vqe_lambda = VQE(ansatz=ansatz, optimizer=optimizer, sampler=sampler)
res_lambda = vqe_lambda.compute_minimum_eigenvalue(operator=L_op)
print("拉普拉斯最小非零特征值 λ1 ≈", res_lambda.eigenvalue.real)

# ==================== 6. QAOA (社区检测) ====================
from qiskit.circuit.library import TwoLocal
from qiskit.algorithms import QAOA
from qiskit.algorithms.optimizers import COBYLA

def modularity_hamiltonian(num_qubits, W):
    k = np.sum(W, axis=1)
    m2 = np.sum(k)
    H = 0
    for i in range(num_qubits):
        for j in range(i+1, num_qubits):
            coeff = -(W[i,j] - k[i]*k[j]/m2)
            if coeff==0: continue
            term = coeff * (Z if i==0 else I^(i-1) ^ Z ^ I^(num_qubits-i-1))
            term = term.compose( Z if j==0 else I^(j-1) ^ Z ^ I^(num_qubits-j-1) )
            H += term
    return PauliSumOp.from_operator(H.to_matrix())

mod_ham = modularity_hamiltonian(num_qubits, W_mat)

p = 2
qaoa_ansatz = TwoLocal(num_qubits, ['ry','rz'], 'cz', entanglement='full', reps=p)
qaoa = QAOA(optimizer=COBYLA(maxiter=200), reps=p,
            quantum_instance=Aer.get_backend('aer_simulator_statevector'),
            ansatz=qaoa_ansatz)

qaoa_res = qaoa.compute_minimum_eigenvalue(operator=mod_ham)
print("Modularity 最小能量 (对应最大社区划分):", qaoa_res.eigenvalue.real)
# 解析二进制解 (0/1 表示社区)
solution = qaoa_res.eigenstate
print("量子社区划分 bitstring:", solution)

# ==================== END ====================
