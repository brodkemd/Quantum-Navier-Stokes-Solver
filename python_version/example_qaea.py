from qiskit.algorithms import EstimationProblem, AmplitudeEstimation

from qiskit.circuit import QuantumCircuit  
from qiskit.utils import QuantumInstance
from qiskit import BasicAer, transpile

import matplotlib.pyplot as plt
import numpy as np


class BernoulliA(QuantumCircuit):
    """A circuit representing the Bernoulli A operator."""

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(theta_p, 0)


class BernoulliQ(QuantumCircuit):
    """A circuit representing the Bernoulli Q operator."""

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        self._theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(2 * self._theta_p, 0)

    def power(self, k):
        # implement the efficient power of Q
        q_k = QuantumCircuit(1)
        print("k =", k , " theta_p=", self._theta_p)
        print("result =", 2 * k * self._theta_p)
        q_k.ry(2 * k * self._theta_p, 0)
        return q_k


p = 0.2

A = BernoulliA(p)
Q = BernoulliQ(p)

problem = EstimationProblem(
    state_preparation=A,  # A operator
    grover_operator=Q,  # Q operator
    objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
)

backend = BasicAer.get_backend("statevector_simulator")
quantum_instance = QuantumInstance(backend)

ae = AmplitudeEstimation(
    num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
    quantum_instance=quantum_instance,
)

ae_result = ae.estimate(problem)
print(ae_result.estimation)
print("Interpolated MLE estimator:", ae_result.mle)

ae_circuit = ae.construct_circuit(problem)
#ae_circuit.decompose().draw(
#    "mpl", style="iqx"
#)

basis_gates = ["h", "ry", "cry", "cx", "ccx", "p", "cp", "x", "s", "sdg", "y", "t", "cz"]
transpile(ae_circuit, basis_gates=basis_gates, optimization_level=2).draw("mpl", style="iqx")

plt.show()