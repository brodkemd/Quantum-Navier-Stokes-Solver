from qiskit import QuantumCircuit, Aer, BasicAer, execute
from qiskit.visualization import(
  plot_state_city,
  plot_bloch_multivector,
  plot_state_paulivec,
  plot_state_hinton,
  plot_state_qsphere)

import matplotlib.pyplot as plt

def show(qc):
    print(qc.draw())
    # execute the quantum circuit
    backend = BasicAer.get_backend('statevector_simulator') # the device to run on
    result = execute(qc, backend).result()
    psi  = result.get_statevector(qc)

    plot_state_qsphere(psi)
    print(psi)
    #state = Statevector.from_instruction(qc)
    #print(state)

simulator = Aer.get_backend('qasm_simulator')

qc = QuantumCircuit(2, 2)
qc.x(1)
qc.h(0)
#qc.cx(1, 0)
#qc.measure(0, 0)
#qc.x(0)
#qc.x(0)
show(qc)
qc.measure(0, 0)
qc.measure(1, 1)

result = execute(qc, backend=simulator, shots=1024).result()
counts  = result.get_counts(qc)
print(counts)
plt.show()

#qc.cx(0, 1)
#show(qc)
