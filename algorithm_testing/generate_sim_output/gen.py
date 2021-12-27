
# qiskit imports
from qiskit.algorithms import AmplitudeEstimation
from qiskit.algorithms import EstimationProblem
from qiskit.circuit import QuantumCircuit
from qiskit import Aer

import numpy as np
import threading
import os
import re


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


class generator:
    # where do you want the code to start in list of data
    start_at_line = 0

    # number of threads this program will utilize
    num_threads = 15

    # counts how many threads are currently done running, init to number of threads
    done_count = num_threads

    # storage lists
    results = []
    thread_input = []
    thread_output = []

    def __init__(self):
        self.ae = AmplitudeEstimation(
            num_eval_qubits=3,  # the number of evaluation qubits specifies circuit width and accuracy
            quantum_instance=Aer.get_backend('qasm_simulator')
        )

        # reading the values from the file
        with open(os.path.join("..", "omega_database"), "r") as f:
            # splits the file's contents at newlines and spaces
            self.nums = re.split('[ \n]', f.read()) 

        # cleans up the data
        self.nums = list(filter(("").__ne__, self.nums))
        
        # converts to a list of floats        
        self.nums = list(map(float, self.nums)) 

        # remembers number of data points
        self.count = len(self.nums)

        input(f"Done loading, going to evaluate {self.count} numbers, press enter to start threads")
        
        # inits threads
        for i in range(self.num_threads):
            # init data communication lists
            self.thread_output.append(0)
            self.thread_input.append(0)

            # creating thread
            print("starting thread", i)
            thread = threading.Thread(target=self.thread_func, args=(i, ))
            thread.daemon = True
            thread.start()

        input(f"Done starting threads, press enter to run")

        # running the main function
        self.Main()

        # recording the results of this code to the proper file
        self.record_result()


    def thread_func(self, id):
        # this function runs on a thread simulates a quantum circuit
        while True:
            # evals to true if main thread reset the numbers in the input list
            if not self.done_count:
                A = BernoulliA(self.thread_input[id])
                Q = BernoulliQ(self.thread_input[id])

                problem = EstimationProblem(
                    state_preparation=A,  # A operator
                    grover_operator=Q,  # Q operator
                    objective_qubits=[0],  # the "good" state Psi1 is identified as measuring |1> in qubit 0
                )
                
                ae_result = self.ae.estimate(problem)
                
                # records this threads result
                self.thread_output[id] = ae_result.mle
                
                # records that this thread has finished
                self.done_count += 1


    def Main(self):
        print("Running...")
        #try:
        for i in range(self.count):
            # setting the value for the thread to run, in order
            self.thread_input[self.num_threads - self.done_count] = self.nums[i]
            
            # indicates that the threads data has been input
            self.done_count -=1

            print(i, "/", self.count, end="\r")

            # if all threads are running
            if not self.done_count:
                # waits for the threads to finish if they have all recieved jobs
                while self.num_threads - self.done_count: pass

                # adds the results of the threads to the results list
                for j in range(len(self.thread_output)): self.results.append(self.thread_output[j])

                    #print("results=", self.results)
        
        # if the user terminated the program
        #except KeyboardInterrupt: pass
        
        print("Recording then Exiting")


    def record_result(self):
        # converts to a list of strings
        self.results = list(map(str, self.results))
        
        # recording in a file
        with open(os.path.join("..", "omega_result_database"), 'w') as f:
            for val in self.results:
                f.write(val + "\n")



inst = generator()