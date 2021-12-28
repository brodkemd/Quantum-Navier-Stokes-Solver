import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import os
import re

class visualize:
    # universal x coordinate system
    x = np.linspace(0, 3, 61)
    vals = []

    def __init__(self):
        # getting the info about the location of this script
        self.head_dir = os.path.split(os.getcwd())[0]
        this_dir = os.path.split(os.getcwd())[1]

        # files that need to be opened, each element gets there own subplot
        files_to_open = [["TempDvals", "TempEvals"], ['MrhoDvals', 'MrhoEvals'], ['MachDvals', 'MachEvals'], ["omega_log"]]
        
        plot_functions = [self.temperature, self.massden, self.machnum, self.none]

        # names of files that indicate how the code was run
        indicating_files = ["QASM", "ORIGINAL"]

        # going through the items in the head direcotry, where the data directories are
        for item in os.listdir(self.head_dir):
            self.item = item # setting so can be used else where
            # if the item is a direcotry and is not the one that this script is in
            if os.path.isdir(os.path.join(self.head_dir, item)) and item != this_dir:
                # looks for the file that indicates how the code was run, if there, name of the file is added to the
                # title of the figure
                indicator = ""
                for indicating_file in indicating_files:
                    if indicating_file in os.listdir(os.path.join(self.head_dir, item)):
                        indicator = indicating_file + "-"

                plt.figure(indicator + item, tight_layout=True, figsize=(15, 8)) # naming the figure

                # counts the subplot
                self.count = 0

                self.num_panes = len(plot_functions) - plot_functions.count(self.none)

                # creates a figure for every element in the in the files_to_open list
                for i, files in enumerate(files_to_open):
                    if plot_functions[i] != self.none:
                        self.count+=1

                    # creates a 1 row many column subplot
                    plt.subplot(2, self.num_panes, self.count)

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

        #self.omega_data_handler()
        plt.show()


    # method for omega_log file
    def omega_log(self, data):
        # adds the values from the inputted list to class list if there isn't one there already
        for data_point in data: self.vals.append(data_point)

        for i, item in enumerate(data):
            data[i] = item.split("->")
            data[i] = list(map(float, data[i]))

        omegas = np.array(data)[:, 0]
        results = np.array(data)[:, 1]
        
        self.count+=1

        # creates a 1 row many column subplot
        plt.subplot(2, self.num_panes, self.count)

        plt.hist(omegas, bins=30)
        plt.title('Distribution of Omega Values')
        plt.xlabel('Omega Value')
        plt.ylabel('number of occurences')

        self.count+=1

        # creates a 1 row many column subplot
        plt.subplot(2, self.num_panes, self.count)

        plt.hist(results, bins=30)
        plt.title('Distribution of Result Values')
        plt.xlabel('Result Value')
        plt.ylabel('number of occurences')

        mean = np.mean(omegas)
        median = np.median(omegas)
        mode = stats.mode(omegas)[0][0]
        std = np.std(omegas)
        variance = np.var(omegas)
        min = np.min(omegas)
        max = np.max(omegas)

        text = f'''Omega:\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}\nstandard deviation = {std:.6f}\nVariance = {variance:.6f}'''

        mean = np.mean(results)
        median = np.median(results)
        mode = stats.mode(results)[0][0]
        std = np.std(results)
        variance = np.var(results)
        min = np.min(results)
        max = np.max(results)

        text+=f'''\n\nResult:\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}\nstandard deviation = {std:.6f}\nVariance = {variance:.6f}'''
        
        plt.figtext(0.7, 0.1, text)

    # for temperature plot
    def temperature(self):
        plt.title('Temperature vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Temperature T')
        plt.legend(['Simulation', "Exact"])

    # method for TempDvals file
    def TempDvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data, "-")
    
    # method for TempEvals file
    def TempEvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)

    # method for mass density plot
    def massden(self):
        plt.title('Mass density vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Mass density mrho')
        plt.legend(['Simulation', "Exact"])

    # method for MrhoDvals file
    def MrhoDvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
    
    # method for MrhoEvals file
    def MrhoEvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)

    # method for mach number plot
    def machnum(self):
        plt.title('Mach number vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Mach number M')
        plt.legend(['Simulation', "Exact"])

    # method for MachDvals file
    def MachDvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
    
    # method for MachEvals file
    def MachEvals(self, data):
        data = list(map(float, data)) # converts to a list of floats
        plt.plot(self.x, data)
    
    def none(self): pass # does nothing, here to make the code general

    def omega_data_handler(self):
        print("Generating Omega data")

        with open(os.path.join(self.head_dir, "database"), 'r') as f:
            nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces
            nums = list(filter(("").__ne__, nums)) # cleans up the data
        
        # adds the values from the inputted list to class list if there isn't one there already
        for num in nums: self.vals.append(num)

        # recording the updated omega values
        with open(os.path.join(self.head_dir, "database"), 'w') as f:
            for val in self.vals:
                f.write(val+"\n")

        for i, item in enumerate(self.vals):
            self.vals[i] = item.split("->")
            self.vals[i] = list(map(float, self.vals[i]))
       
        omegas = np.array(self.vals)[:, 0]
        results = np.array(self.vals)[:, 1]

        # figure that displays information on the omega database
        plt.figure("Omega database", tight_layout=True, figsize=(12, 5))

        plt.subplot(1, 3, 1)
        plt.hist(omegas, bins=30)
        plt.title('Distribution of Omega Values')
        plt.xlabel('Omega Value')
        plt.ylabel('number of occurences')
        
        plt.subplot(1, 3, 2)
        plt.hist(results, bins=30)
        plt.title('Distribution of Result Values')
        plt.xlabel('Result Value')
        plt.ylabel('number of occurences')
        
        '''
        x = list(np.linspace(0, 1, len(self.vals)))
        plt.plot(x, omegas, '.')
        plt.plot(x, results, '.')
        plt.title('Value Progression Through Time')
        plt.ylabel('Value')
        plt.xlabel('Position of occurence in the progression')
        plt.legend(['Omega', "Estimate"])
        '''

        mean = np.mean(omegas)
        median = np.median(omegas)
        mode = stats.mode(omegas)[0][0]
        std = np.std(omegas)
        variance = np.var(omegas)
        min = np.min(omegas)
        max = np.max(omegas)

        text = f'''Omega:\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}\nstandard deviation = {std:.6f}\nVariance = {variance:.6f}'''

        mean = np.mean(results)
        median = np.median(results)
        mode = stats.mode(results)[0][0]
        std = np.std(results)
        variance = np.var(results)
        min = np.min(results)
        max = np.max(results)

        text+=f'''\nResult:\nMean = {mean:.6f}\nMedian = {median:.6f}\nMode = {mode:.6f}\nMin = {min:.6f}\nMax = {max:.6f}\nRange = {(max - min):.6f}\nstandard deviation = {std:.6f}\nVariance = {variance:.6f}'''

        plt.figtext(0.7, 0.1, text)
        
visualize()
