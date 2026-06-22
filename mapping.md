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
