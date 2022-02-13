from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit, Aer
import matplotlib.pyplot as plt
import numpy as np
import os

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
        self.strings.append("qc." + string)

    def __init__(self, p, num_eval_qbits, num_control_qbits):
        self.num_eval_qbits = num_eval_qbits + 1
        self.num_control_qbits = num_control_qbits
        print(f"{num_eval_qbits + num_control_qbits}")
        self.strings = [f"from qiskit import QuantumCircuit", 
                        f"qc = QuantumCircuit({num_eval_qbits + num_control_qbits})"]

        # inits the qubits
        for i in range(1, self.num_eval_qbits):
            self.add(f"h({i})")            

        theta_p = 2 * np.arcsin(np.sqrt(p))
        self.add(f"ry({theta_p}, 0)")
        if self.num_control_qbits > 1:
            self.add(f"ry({theta_p}, {self.num_control_qbits + self.num_eval_qbits - 1})")

        # running the function that generates the main portion of the algorithm
        self.loop()


    def repeat_gates(self):
        '''
        creates gates that are repeatedly used in the algorithm

        '''
        self.add("sdg(0)")
        self.add("h(0)")
        self.add("sdg(0)")

    
    def end_gates(self):
        '''
        gates that finsish the algorithm
        
        '''
        for i in range(1, self.num_eval_qbits - 1):
            self.add(f"h({i})")
            for j in range(i + 1, self.num_eval_qbits):
                self.add(f"cp({-1 * np.pi / (2**(j - i))}, {i}, {j})")
            
        self.add(f"h({self.num_eval_qbits - 1})")


    def loop(self):
        '''
        creates the majority of the gates
        
        '''
        count = self.num_eval_qbits - 1
        b = False

        vals = [2.21429743558818, 4.06888787159, 1.28700221758657, 4.99618308959302, -0.567588218416656, 6.85077352559624, -4.27676909042310, 10.5599543976027]
        #vals.reverse()
        itr = 0
        while True:
            self.add(f"cx({int(np.ceil(count))}, 0)")

            self.repeat_gates()
            self.add(f"p({vals[itr]}, 0)")
            
            self.repeat_gates()
            self.add(f"p({3 * np.pi}, 0)")
                        
            if count < 1: break
            count-=0.5
            b = not(b)
            itr+=1
        
        # adds the gates to the end
        self.end_gates()


    def compile_to_file(self, filename, include_statevector=False, show_after=False):
        print("Compiling to:", filename)
        with open(filename, 'w') as f:
            if show_after:
                self.strings.insert(0, "import matplotlib.pyplot as plt")
            for string in self.strings:
                f.write(string + "\n")
                if include_statevector:
                    exec(string)
                    if "qc." in string:
                        exec("result = simulated_state_vector(qc)")
                        exec("list_string = \"#Operation resulted in the statevector:\\n#[\"")
                        exec("for item in result.get_statevector(qc, decimals=3): list_string+=(str(item)+\" \")")
                        exec("f.write(list_string + \"]\\n\")")

                f.write("\n")
            
            if show_after:
                f.write("qc.draw(output=\"mpl\")\nplt.tight_layout()\nplt.show()")
                os.system(f"python {filename}")


    def disp(self):
        '''
        draws the circuit in a matplotlib window then displays it

        '''
        self.draw(output='mpl')
        plt.show()

p = 0.2
num_eval_qbits = 4
num_control_qbits = 1

# running the above class
inst = qaea(p, num_eval_qbits, 1)
inst.compile_to_file("out1.py", include_statevector=True, show_after=True)