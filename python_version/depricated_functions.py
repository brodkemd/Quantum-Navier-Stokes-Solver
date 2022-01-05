from pprint import pprint
import numpy as np
import os

def QAmpEst(self, omega):
        Estimates = np.zeros((int(self.TotRuns)))

        # start loop to carry out TotRuns simulation runs
        for runs in range(int(self.TotRuns)):
            randev = self.randQAEA(omega)    # randQAEA generates random deviate 
                                        # with probability distribution
                                        # produced by QAEA
            Estimates[runs] = randev

        return (np.sin(np.pi*np.median(Estimates)/self.N))**(2)

def MeanOrc(self, Gij):
    #MEANORC is an oracle function for mean value of g_ij.
    temp = 0      # used to accumulate mean value

    #  accumulate mean value
    for j in range(int(self.N)): temp = temp + Gij[j]

    return temp/self.N    # RHS is mean value


def randQAEA(self, omega):
    #randQAEA Generates random deviate for Quantum Amplitude Estimation.

    Momega = self.N*omega

    Tiny = 1*10**(-1 * 50)       # tiny number used to prevent 0/0 in pofx calculation

    x = -1
    ratio = -1

    # begin calculation of randev
    while ((x < 0 or x > (self.N-1)) or (np.random.rand(1) > ratio)):
        v1 = np.random.rand(1)
        v2 = 2*np.random.rand(1) - 1
        magv = v1**2 + v2**2

        while (magv > 1):
            v1 = np.random.rand(1)
            v2 = 2*np.random.rand(1) - 1
            magv = v1**2 + v2**2

        y = v2/v1
        x = y + Momega
        intx = np.round(x)
        inty = intx - Momega

        if (intx >= 0):
            if (intx <= (self.N-1)):
                nearx = intx
            else:
                nearx = self.N

        else:
            # 'Warning: intx is negative - set nearx to intx'
            nearx = intx


        # pofx evaluates QAEA probability distribution at nearx
        pofx = (1/2)*(np.sin(np.pi*(Momega - nearx + Tiny)))**(2)/(self.N*np.sin((np.pi/self.N)*(Momega - nearx + Tiny)))**(2)+(1/2)*(np.sin(np.pi*(self.N - Momega - nearx + Tiny)))**(2)/(self.N*np.sin((np.pi/self.N)*(self.N - Momega - nearx + Tiny)))**(2)
        ratio = (1 + (inty)**(2))*pofx

    # nearx is the desired random deviate randev
    return nearx

def FuncOrc(self, t, TCoeffs, rmaxp1):
    #FUNCORC evaluates g_ij[u] at N knot times in subsubinterval j
    # initialize parameters and arrays
    # f stores d components of ODE driver function f at each interior 
    #    grid-point and N knot times for subsubinterval j
    f = np.zeros((int(self.d), int(self.Tot_Int_Pts), int(self.N)))   

    # evaluate f at N knot times for subsubinterval j and each interior
    # grid-point
    #self.stop("VAl", self.fOrc(t[3], t[0], TCoeffs, rmaxp1))
    for k in range(int(self.N)): f[:, :, k] = self.fOrc(t[k], t[0], TCoeffs, rmaxp1)

    # assign values of Gij at N knot times t for subsubinterval j and each 
    #       interior grid-point
    return f

def stop(self, name, val):
    print(name)
    pprint(val)
    exit(0)


def save(self):
    with open("self_vals", 'w') as f:
        for item in sorted(dir(self), key=str.lower):
            if "__" not in item and not callable(eval("self." +item)):
                if isinstance(eval("self." + item), (np.ndarray, list)):
                    f.write(item + " =\n")
                    np.savetxt(f, eval("self." + item), "%f", delimiter=" ")
                else:
                    f.write(item + " = " + str(eval("self." + item)) + "\n")
                f.write("\n")


def check(self, items_in_list):
    check_l = []
    for item in dir(self):
        if "__" not in item and not callable(eval("self." +item)):
            check_l.append(item)
    
    items_in_list = items_in_list.split(", ")

    for item in items_in_list:
        if item not in check_l:
            raise ValueError(item + " not in variables")


def show(self, val=True):
    for item in sorted(dir(self), key=str.lower):
        if "__" not in item and not callable(eval("self." +item)):
            if val: print(item, "=", eval("self." + item))
            else: print(item)

def cleanup():
    for item in os.listdir('.'):
        if os.path.isfile(item) and '.' not in item:
            print(item)
            option = int(input("Good to remove (0=no, 1=yes):"))
            if option:
                print("removing")
                os.remove(item)
            print()
    
    exit()