import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import os, re

# font control
NORMAL_SIZE = 17
BIGGER_SIZE = 20

plt.rc('font', size=NORMAL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=NORMAL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=NORMAL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=NORMAL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=NORMAL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=NORMAL_SIZE)    # legend fontsize

class visualize:
    # universal x coordinate system
    x = np.linspace(0, 3, 61)
    vals = []
    differnce_vals = []
    legend = []

    def __init__(self):
        # getting the info about the location of this script
        self.head_dir = os.path.split(os.getcwd())[0]

        ignore_these = [os.path.split(os.getcwd())[1], "c++_backend", "orig_matlab_files", "orig_python_files", "test", "__pycache__"]

        # files that need to be opened, each element gets there own subplot
        files_to_open = [["Temp_D", "Temp_E"], ['Mrho_D', 'Mrho_E'], ['Mach_D', 'Mach_E'], ["Press_D", "Press_E"]]
        
        plot_functions = [self.temperature, self.massden, self.machnum, self.pressnum]

        # names of files that indicate how the code was run
        indicating_files = ["QASM", "ORIGINAL", "REAL", "INTERPOLATED", "ESTIMATED"]

        # going through the items in the head direcotry, where the data directories are
        for item in os.listdir(self.head_dir):
            self.item = item # setting so can be used else where
            # if the item is a direcotry and is not the one that this script is in
            if os.path.isdir(os.path.join(self.head_dir, item)) and item not in ignore_these:
                # looks for the file that indicates how the code was run, if there, name of the file is added to the
                # title of the figure
                indicator = ""
                for indicating_file in indicating_files:
                    if indicating_file in os.listdir(os.path.join(self.head_dir, item)):
                        indicator+=(indicating_file + "-")

                plt.figure(indicator + item, tight_layout=True, figsize=(15, 8)) # naming the figure

                # counts the subplot
                self.count = 0

                self.num_panes = 2#len(plot_functions) - plot_functions.count(self.none)

                # creates a figure for every element in the in the files_to_open list
                for i, files in enumerate(files_to_open):
                    if plot_functions[i] != self.none:
                        self.count+=1

                    # creates a 1 row many column subplot
                    plt.subplot(self.num_panes, 2, self.count)

                    # plots all of the file names in the list at the current position in files_to_open on the same figure
                    for file in files:
                        try:
                            # opens the file
                            print("opening", os.path.join(self.head_dir, item, file))
                            with open(os.path.join(self.head_dir, item, file), 'r') as f:
                                nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces
                                nums = list(filter(("").__ne__, nums)) # cleans up the data
                                                                
                                # runs the functions plotting function
                                exec(f"self.{file}({nums})")
                        
                        # catches if the file does not exist, tells user
                        except FileNotFoundError as e: print(e)

                    # runs the files formating
                    plot_functions[i]()
                
                plt.savefig((indicator + item).replace(":", ""))
        #plt.show()



    # for temperature plot
    def temperature(self):
        plt.title('Temperature vs. nozzle position', fontsize = BIGGER_SIZE)
        plt.xlabel('Nozzle position')
        plt.ylabel('Temperature')
        #plt.legend(['Simulation', "1-D analytic solution"])
        self.legend.append("Temperature")

    # method for TempDvals file
    def Temp_D(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data, "-")
        self.differnce_vals.append(data)
    
    # method for TempEvals file
    def Temp_E(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)

    # method for mass density plot
    def massden(self):
        plt.title('Mass density vs. nozzle position', fontsize = BIGGER_SIZE)
        plt.xlabel('Nozzle position')
        plt.ylabel('Mass density')
        #plt.legend(['Simulation', "1-D analytic solution"])
        self.legend.append("Mass Density")

    # method for MrhoDvals file
    def Mrho_D(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)
    
    # method for MrhoEvals file
    def Mrho_E(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)

    # method for mach number plot
    def machnum(self):
        plt.title('Mach number vs. nozzle position', fontsize = BIGGER_SIZE)
        plt.xlabel('Nozzle position')
        plt.ylabel('Mach number')
        #plt.legend(['Simulation', "1-D analytic solution"])
        self.legend.append("Mach Number")

    # method for MachDvals file
    def Mach_D(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)
    
    # method for MachEvals file
    def Mach_E(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)

    def pressnum(self):
        plt.title('Pressure vs. nozzle position', fontsize = BIGGER_SIZE)
        plt.xlabel('Nozzle position')
        plt.ylabel('Pressure')
        #plt.legend(['Simulation', "1-D analytic solution"])
        self.legend.append("Pressure")

    # method for MachDvals file
    def Press_D(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)

    # method for MachEvals file
    def Press_E(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
        self.differnce_vals.append(data)
    
    def none(self): pass # does nothing, here to make the code general

visualize()


"""
    def error(self):
        self.count+=1

        # creates a 1 row many column subplot
        plt.subplot(self.num_panes, 2, self.count)

        plt.title('Error in simulated values')
        plt.xlabel('nozzle position')
        plt.ylabel('error in value')

        for i in range(0, len(self.differnce_vals), 2):
            plt.plot(self.x, np.subtract(self.differnce_vals[i], self.differnce_vals[i+1]))

        plt.legend(self.legend)

        self.legend.clear()
        self.differnce_vals.clear()

    # method for omega_log file
    def QAmpEst(self, data):
        return
        # adds the values from the inputted list to class list if there isn't one there already
        for data_point in data: self.vals.append(data_point)

        for i, item in enumerate(data):
            data[i] = item.split("->")
            data[i] = list(map(float, data[i]))

        omegas = np.array(data)[:, 0]
        results = np.array(data)[:, 1]
        
        self.error()

        self.count+=1

        # creates a 1 row many column subplot
        plt.subplot(self.num_panes, 2, self.count)

        plt.hist(omegas, bins=30)
        plt.title('Distribution of Omega Values')
        plt.xlabel('Omega Value')
        plt.ylabel('number of occurences')

        self.count+=1

        # creates a 1 row many column subplot
        plt.subplot(self.num_panes, 2, self.count)

        plt.hist(results, bins=30)
        plt.title('Distribution of Result Values')
        plt.xlabel('Result Value')
        plt.ylabel('number of occurences')
        
        return
        
        mean = np.mean(omegas)
        median = np.median(omegas)
        mode = stats.mode(omegas)[0][0]

        min = np.min(omegas)
        max = np.max(omegas)

        text = f'''Omega:\nCount={len(omegas)}\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}'''

        mean = np.mean(results)
        median = np.median(results)
        mode = stats.mode(results)[0][0]
        min = np.min(results)
        max = np.max(results)

        text+=f'''\n\nResult:\nCount={len(results)}\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}'''
        
        plt.figtext(0.8, 0.1, text)
"""