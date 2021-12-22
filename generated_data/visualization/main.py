import matplotlib.pyplot as plt
import numpy as np
import os

class visualize:
    # universal x coordinate system
    x = np.linspace(0, 3, 61)


    def __init__(self):
        # getting the info about the location of this script
        head_dir = os.path.split(os.getcwd())[0]
        this_dir = os.path.split(os.getcwd())[1]

        # files that need to be opened
        files_to_open = [["TempDvals", "TempEvals"], ['MrhoDvals', 'MrhoEvals'], ['MachDvals', 'MachEvals']]
        
        plot_functions = [self.temperature, self.massden, self.machnum]

        indicating_files = ["QASM", "ORIGINAL"]

        # going through the items in the head direcotry, where the data directories are
        for item in os.listdir(head_dir):
            # if the item is a direcotry and is not the one that this script is in
            if os.path.isdir(os.path.join(head_dir, item)) and item != this_dir:
                indicator = ""
                for indicating_file in indicating_files:
                    if indicating_file in os.listdir(os.path.join(head_dir, item)):
                        indicator = indicating_file + "-"

                plt.figure(indicator + item) # naming the figure

                count = 0

                # creates a figure for every element in the in the files_to_open list
                for i, files in enumerate(files_to_open):
                    count+=1

                    plt.subplot(1, len(plot_functions), count)

                    # plots all of the file names in the list at the current position in files_to_open on the same figure
                    for file in files:
                        try:
                            # opens the file
                            with open(os.path.join(head_dir, item, file), 'r') as f:
                                # gets the data from the file in a useable format
                                nums = f.read().split()
                                nums = list(map(float, nums))
                                
                                # runs the functions plotting function in file_methods.py
                                exec(f"self.{file}({nums})")
                        
                        # catches if the file does not exist, tells user
                        except FileNotFoundError:
                            print(f"ERROR: Unable to open {file} from {os.path.join(head_dir, item)}")

                    # runs the files formating
                    plot_functions[i]()

        plt.show()

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

visualize()
