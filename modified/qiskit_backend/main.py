import time
start = time.time()

import socket
import struct

# qiskit imports
from qiskit.algorithms import AmplitudeEstimation
from qiskit.algorithms import EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit.utils import QuantumInstance
from qiskit import transpile
from qiskit import BasicAer
from qiskit import Aer

import numpy as np

# create a socket object
serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

port = 9995

# bind to the port
serversocket.bind(('localhost', port))

# queue 1 request
serversocket.listen(1)

print("initialized server with port number:", port)

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


while True:
    # establish a connection
    clientsocket, addr = serversocket.accept()

    #print("Connected to %s" % str(addr))

    data = clientsocket.recv(4) # recieves a packed bytes
    p = struct.unpack('f', data)[0]
    #print("recieved:", p)
    
    A = BernoulliA(p)
    Q = BernoulliQ(p)

    problem = EstimationProblem(
        state_preparation=A,  # A operator
        grover_operator=Q,  # Q operator
        objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
    )

    ae_result = ae.estimate(problem)
    #print('sending:', ae_result.mle)

    data = struct.pack('f', ae_result.mle)

    clientsocket.sendall(data)
    clientsocket.close()