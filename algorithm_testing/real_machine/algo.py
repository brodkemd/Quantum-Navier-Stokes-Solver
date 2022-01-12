nums = []

# qiskit imports
from qiskit.algorithms import AmplitudeEstimation, EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit import IBMQ, Aer
import numpy as np
import time, warnings

warnings.filterwarnings("ignore")

# loading my account
IBMQ.load_account()

# telling it where to look for the quantum devices
backend = IBMQ.get_provider('ibm-q').get_backend('ibmq_belem')


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


ae = AmplitudeEstimation(
    num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
    quantum_instance=backend
)

omega = 0.25
A = BernoulliA(omega)
Q = BernoulliQ(omega)

problem = EstimationProblem(
state_preparation=A,  # A operator
grover_operator=Q,  # Q operator
objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
)

ae = AmplitudeEstimation(
    num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
    quantum_instance=backend
)

print("running")
result = ae.estimate(problem)
print("estimate=", result.estimation)
print("interpolate=", result.mle)

exit()
# iteration counter
i = 0
# loops through all of the devices available to me
for device in provider.backends():
    # this excludes simulators from the list of real devices, this is good
    if device.name().count("simulator") != 0 or device.name() == "ibmq_armonk": continue

    # adding the current device to the list
    devices[device] = device.status().pending_jobs

    # printing the iteration of the loop to be used as an option for the user later, the device name, 
    # how many other jobs are in the queue, and the number of Qbits each device can handle
    print(f"{i} : {device.name()} has {device.status().pending_jobs} queued")
    print("Qbits=", len(device.properties().qubits))
    print("Gates=", device.properties().gates)

    i = i + 1

# essentially asking the user what device they want
#option = int(input("What device? (input the integer before the \":\" in each line): "))
option = min(devices, key=devices.get)

# getting the device the user wanted
print("using:", option)

exit()

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

print("running", end=" ")

if not real_machine:
    print("on qasm simulator")
    log_file = f"qasm_{time.ctime().replace(' ', '_')}"
    backend = Aer.get_backend('qasm_simulator')

    ae = AmplitudeEstimation(
        num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
        quantum_instance=backend
    )

else:
    print("on real machine")
    log_file = f"real_{time.ctime().replace(' ', '_')}"
    # loading my account
    IBMQ.load_account()

    # telling it where to look for the quantum devices
    provider = IBMQ.get_provider('ibm-q')

for omega in nums:
    A = BernoulliA(omega)
    Q = BernoulliQ(omega)

    problem = EstimationProblem(
        state_preparation=A,  # A operator
        grover_operator=Q,  # Q operator
        objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
    )

    if real_machine:
        # list that holds all of the names of the devices
        devices = {}

        # iteration counter
        i = 0
        # loops through all of the devices available to me
        for device in provider.backends():
            # this excludes simulators from the list of real devices, this is good
            if device.name().count("simulator") != 0 or device.name() == "ibmq_armonk": continue

            # adding the current device to the list
            devices[device] = device.status().pending_jobs

            # printing the iteration of the loop to be used as an option for the user later, the device name, 
            # how many other jobs are in the queue, and the number of Qbits each device can handle
            #print(f"{i} : {device.name()} has {device.status().pending_jobs} queued")

            #i = i + 1

        # essentially asking the user what device they want
        #option = int(input("What device? (input the integer before the \":\" in each line): "))
        option = min(devices, key=devices.get)

        # getting the device the user wanted
        print("using:", option)

        #input("good?")

        backend = provider.get_backend(option.name().strip()) # for some reason the name of the device has extra whitespaces so I strip them

        ae = AmplitudeEstimation(
            num_eval_qubits=4,  # the number of evaluation qubits specifies circuit width and accuracy
            quantum_instance=backend
        )

    result = (np.sin(np.pi * ae.estimate(problem).mle))**2

    # recording the result
    with open(log_file, "a") as f: f.write(str(result) + "\n")