import time
start = time.time()

# qiskit imports
from qiskit.algorithms import AmplitudeEstimation
from qiskit.algorithms import EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit import Aer
import numpy as np

# data communication
import socket
import struct
import os

# logging functions
import logger

# if the values sent to this code should be recorded in a log file
log = True

# if the result should be computed or not
compute = False

# true if you want to record to the result to the omega value

# names the file to be made based on the time
log_file = f"omega_log"

# create a socket object
serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
#serversocket.setblocking(False)

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
    num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
    quantum_instance=backend
)

print("load time:", time.time() - start)

# init val that shows when the server should disconnect from the client
close = False

while True:
    # establish a connection
    clientsocket, addr = serversocket.accept()

    print("Connected to %s" % str(addr))
  
    # creates/erases the log file
    with open(log_file, "w") as f: pass

    # receives two config bytes from the client
    config = clientsocket.recv(2)
    print(f"len = {len(config)}")
    options = struct.unpack("??", config)
    
    # parsing the options to tell the server operator their meaning
    if options[0]: opt0 = "run the simulator"
    else: opt0 = "not run anything"
    if options[1]: opt1 = "log"
    else: opt1 = "not log anything"
    print(f"Will {opt0} and {opt1}")
    input("press enter to continue")

    # reads forever until a timeout
    while True:
        # if the server should stop, no data to log at this point
        if close:
            logger.Main(log_file, simulated=options[0]) # recording the data output
            os.remove(log_file) # removing the log file
            break # breaks inner loop and closes tcp connection

               
        # if should compute the result, then computing and sending back
        if options[0]:
            # recieves a packed bytes
            data = clientsocket.recv(8)

            # if invalid data is returned then it is most likely the end of the program
            try:
                omega = struct.unpack('d', data)[0]
            except struct.error: # will throw this error
                close = True

            A = BernoulliA(omega)
            Q = BernoulliQ(omega)

            problem = EstimationProblem(
                state_preparation=A,  # A operator
                grover_operator=Q,  # Q operator
                objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
            )

            result = ae.estimate(problem).mle
            #print('sending:', ae_result.mle)

            data = struct.pack('d', result)

            clientsocket.sendall(data)
        
        # the client computed then, recording their data
        else:
            # recieves a packed bytes
            data = clientsocket.recv(16)
            data = struct.unpack("dd", data)
            omega = data[0]
            result = data[1]

        if options[1]:
            # if should log, then logging to the log file
            with open(log_file, 'a') as f:
                print('logging:', f"{omega}->{result}\n")
                f.write(f"{omega}->{result}\n")
        
    clientsocket.close()