p = 0.2
num_qbits = 5

if 0:
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
    ae_circuit.decompose().draw(
        "mpl", style="iqx"
    )

    basis_gates = ["h", "ry", "cry", "cx", "ccx", "p", "cp", "x", "s", "sdg", "y", "t", "cz"]
    transpile(ae_circuit, basis_gates=basis_gates, optimization_level=2).draw("mpl", style="iqx")

    plt.show()

else:
    from qiskit.quantum_info import Statevector
    from qiskit import QuantumCircuit, Aer
    import matplotlib.pyplot as plt
    import numpy as np

    def simulated_state_vector(quantum_circuit):
        backend = Aer.get_backend('statevector_simulator')
        job = backend.run(quantum_circuit)
        result = job.result()

        return result

    class qaea(QuantumCircuit):
        '''
        my rendition of quantum amplitude estimation algorithm built into
        qiskit

        should be able to intertwine it better with the rest of ns_q.py and 
        possibly make some upgrades in the future

        '''
        def add(self, string):
            self.strings.append(string.replace("self", "qc"))

        def __init__(self, p, num_qbits):
            self.num_qbits = num_qbits

            self.strings = [f"qc = QuantumCircuit({num_qbits})"]

            # creates the circuit
            super().__init__(self.num_qbits)

            # inits the qubits
            for i in range(1, self.num_qbits):
                self.h(i)
                self.add(f"self.h({i})")
                self.p(0, i)
                self.add(f"self.p(0, {i})")
                

            theta_p = 2 * np.arcsin(np.sqrt(p))
            self.ry(theta_p, 0)
            self.add(f"self.ry({theta_p}, 0)")

            # for debugging
            print("theta_p:", theta_p)
            self.barrier()

            # running the function that generates the main portion of the algorithm
            self.loop()

            self.strings.append("print(qc.draw())")


        def repeat_gates(self):
            '''
            creates gates that are repeatedly used in the algorithm

            '''
            self.sdg(0)
            self.add("self.sdg(0)")
            
            self.h(0)
            self.add("self.h(0)")

            self.sdg(0)
            self.add("self.sdg(0)")
        
        
        def end_gates(self):
            '''
            gates that finsish the algorithm
            
            '''
            for i in range(1, self.num_qbits - 1):
                self.h(i)
                self.add(f"self.h({i})")
                for j in range(i + 1, self.num_qbits):
                    self.cp(-1 * np.pi / (2**(j - i)), i, j)
                    self.add(f"self.cp({-1 * np.pi / (2**(j - i))}, {i}, {j})")
                
                self.barrier()
            
            self.h(self.num_qbits - 1)
            self.add(f"self.h({self.num_qbits - 1})")


        def loop(self):
            '''
            creates the majority of the gates
            
            '''
            count = self.num_qbits - 1
            b = False

            vals = [2.21429743558818, 4.06888787159, 1.28700221758657, 4.99618308959302, -0.567588218416656, 6.85077352559624, -4.27676909042310, 10.5599543976027]

            itr = 0
            while True:
                self.cx(int(np.ceil(count)), 0)
                self.add(f"self.cx({int(np.ceil(count))}, 0)")

                if (b):
                    self.p(0, 0)
                    self.add(f"self.p(0, 0)")

                self.repeat_gates()
                self.p(vals[itr], 0)
                self.add(f"self.p({vals[itr]}, 0)")
                self.repeat_gates()
                self.p(3 * np.pi, 0)
                self.add(f"self.p({3 * np.pi}, 0)")
                
                self.barrier() # for looks
                if count < 1: break
                count-=0.5
                b = not(b)
                itr+=1
            
            # adds the gates to the end
            self.end_gates()


        def run_string(self):
            for item in self.strings:
                print(f"Running: {item}")
                exec(item)

        def show_state_vector_evolution(self):
            import sys
            original_stdout = sys.stdout
            with open('out.txt', 'w') as f:
                sys.stdout = f # Change the standard output to the file we created.
                #print('This message will be written to a file.')
                #sys.stdout = original_stdout 
                for item in self.strings:
                    print(item)
                    exec(item)
                    exec("result = simulated_state_vector(qc)")
                    exec("print(result.get_statevector(qc, decimals=3))")
                    input()

            sys.stdout = original_stdout

        def disp(self):
            '''
            draws the circuit in a matplotlib window then displays it

            '''
            for item in self.strings:
                print(item)
            #self.draw(output='mpl')
            #plt.show()


    # running the above class
    inst = qaea(p, num_qbits)

    # showing the circuit it generates
    #inst.disp()
    #inst.run_string()
    inst.show_state_vector_evolution()