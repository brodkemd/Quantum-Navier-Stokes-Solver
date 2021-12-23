import matplotlib.pyplot as plt
import numpy as np
import os
import re

class visualize:
    # universal x coordinate system
    x = np.linspace(0, 3, 61)
    omega_vals = []

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

                plt.figure(indicator + item, tight_layout=True, figsize=(15, 5)) # naming the figure

                # counts the subplot
                count = 0

                num_panes = len(plot_functions) - plot_functions.count(self.none)

                # creates a figure for every element in the in the files_to_open list
                for i, files in enumerate(files_to_open):
                    if plot_functions[i] != self.none:
                        count+=1

                    # creates a 1 row many column subplot
                    plt.subplot(1, num_panes, count)

                    # plots all of the file names in the list at the current position in files_to_open on the same figure
                    for file in files:
                        try:
                            # opens the file
                            print("opening", os.path.join(self.head_dir, item, file))
                            with open(os.path.join(self.head_dir, item, file), 'r') as f:
                                nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces
                                nums = list(filter(("").__ne__, nums)) # cleans up the data
                                nums = list(map(float, nums)) # converts to a list of floats
                                
                                # runs the functions plotting function
                                exec(f"self.{file}({nums})")
                        
                        # catches if the file does not exist, tells user
                        except FileNotFoundError as e: print(e)

                    # runs the files formating
                    plot_functions[i]()

        self.omega_data_handler()

        plt.show()


    # method for omega_log file
    def omega_log(self, data):
        data = list(set(data)) # removes duplicates
        
        # adds the values from the inputted list to class list if there isn't one there already
        for data_point in data:
            self.omega_vals.append(data_point)

    # for temperature plot
    def temperature(self):
        plt.title('Temperature vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Temperature T')
        plt.legend(['Simulation', "Exact"])

    # method for TempDvals file
    def TempDvals(self, data):
        plt.plot(self.x, data, "-")
    
    # method for TempEvals file
    def TempEvals(self, data):
        plt.plot(self.x, data)

    # method for mass density plot
    def massden(self):
        plt.title('Mass density vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Mass density mrho')
        plt.legend(['Simulation', "Exact"])

    # method for MrhoDvals file
    def MrhoDvals(self, data):
        plt.plot(self.x, data)
    
    # method for MrhoEvals file
    def MrhoEvals(self, data):
        plt.plot(self.x, data)

    # method for mach number plot
    def machnum(self):
        plt.title('Mach number vs. nozzle position')
        plt.xlabel('Nozzle position x')
        plt.ylabel('Mach number M')
        plt.legend(['Simulation', "Exact"])

    # method for MachDvals file
    def MachDvals(self, data):
        plt.plot(self.x, data)
    
    # method for MachEvals file
    def MachEvals(self, data):
        plt.plot(self.x, data)
    
    def none(self): pass # does nothing, here to make the code general

    def omega_data_handler(self):
        print("Generating Omega data")

        with open(os.path.join(self.head_dir, "omega_database"), 'r') as f:
            nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces
            nums = list(filter(("").__ne__, nums)) # cleans up the data
            nums = list(map(float, nums)) # converts to a list of floats
       
        # adds the values from the inputted list to class list if there isn't one there already
        for num in nums: self.omega_vals.append(num)

        # removes duplicates
        self.omega_vals = list(set(self.omega_vals))

        # figure that displays information on the omega database
        plt.figure("Omega database", tight_layout=True, figsize=(12, 5))

        plt.subplot(1, 3, 1)
        plt.hist(self.omega_vals, bins=30)
        plt.title('Distribution of Omega Values')
        plt.xlabel('Omega Value')
        plt.ylabel('number of occurences')
        
        plt.subplot(1, 3, 2)
        x = list(np.linspace(0, 1, len(self.omega_vals)))
        plt.plot(x, self.omega_vals, '.')
        plt.title('Omega Values Progression')
        plt.ylabel('Omega Value')
        plt.xlabel('Position of occurence in the progression')

        mean = np.mean(self.omega_vals)
        median = np.median(self.omega_vals)
        std = np.std(self.omega_vals)
        variance = np.var(self.omega_vals)

        plt.figtext(0.7, 0.7,
f'''
Mean = {mean}\n
Median = {median}\n
standard deviation = {std}\n
Variance = {variance}
''')
        
        self.omega_vals = list(map(str, self.omega_vals)) # converts to a list of strings

        # recording the updated omega values
        with open(os.path.join(self.head_dir, "omega_database"), 'w') as f:
            for omega in self.omega_vals:
                f.write(omega+"\n")

visualize()
