
import time
start = time.time()

# qiskit imports
from qiskit.algorithms import AmplitudeEstimation
from qiskit.algorithms import EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit.utils import QuantumInstance
from qiskit import transpile
from qiskit import BasicAer
from qiskit import Aer

import numpy as np
import sys, os

class BernoulliA(QuantumCircuit):
    """
    A circuit representing the Bernoulli A operator.
    
    """

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(theta_p, 0)


class BernoulliQ(QuantumCircuit):
    """
    A circuit representing the Bernoulli Q operator.
    
    """

    def __init__(self, probability):
        super().__init__(1)  # circuit on 1 qubit

        self._theta_p = 2 * np.arcsin(np.sqrt(probability))
        self.ry(2 * self._theta_p, 0)

    def power(self, k):
        # implement the efficient power of Q
        q_k = QuantumCircuit(1)
        q_k.ry(2 * k * self._theta_p, 0)
        return q_k

backend = Aer.get_backend('qasm_simulator')

ae = AmplitudeEstimation(
    num_eval_qubits=3,  # the number of evaluation qubits specifies circuit width and accuracy
    quantum_instance=backend
)

print("load time:", time.time() - start)

while "STOP" not in os.listdir():
    if "START" not in os.listdir():
        continue
    start = time.time()

    os.remove("START")

    p = get_p()

    A = BernoulliA(p)
    Q = BernoulliQ(p)

    problem = EstimationProblem(
        state_preparation=A,  # A operator
        grover_operator=Q,  # Q operator
        objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
    )

    ae_result = ae.estimate(problem)

    return_estimate(ae_result.mle)

    f = open("DONE", 'w')
    f.close()


    #print(ae_result.estimation)
    #print(ae_result.mle) # interpolated result
    #print("Exe time:", time.time() - start)

os.remove("STOP")
