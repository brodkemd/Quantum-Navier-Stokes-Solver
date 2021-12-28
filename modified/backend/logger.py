import os, time, shutil

def Main(log_file=None, simulated=False):
    # where all data holding directories are stored
    data_storage = os.path.join("..", "..", "generated_data")
    data_origin_directory = ".."
    data_destination = os.path.join(data_storage, f"{time.ctime().replace(' ', '_')}")

    # creates a directory in the correct directory named with the data and time that holds the data on the run
    os.mkdir(data_destination)

    # creating a dict of files to move to log directory, key is the name and value at the key is the relative path
    to_move = {log_file : log_file}
    
    # adds the names of the data files from the directory with them to the list
    for item in os.listdir(data_origin_directory):
        #print(not os.path.isdir(data_origin_directory + item), )
        if (not os.path.isdir(os.path.join(data_origin_directory, item)) and "." not in item) or item == "FINISHED.txt":
            #print(os.path.join(data_origin_directory, item), "->", not os.path.isdir(os.path.join(data_origin_directory, item)))
            to_move[item] = os.path.join(data_origin_directory, item)
    
    # copies the data files to the data directory
    for key in to_move: shutil.copyfile(to_move[key], os.path.join(data_destination, key))
    
    # if the server ran the simulator, creates a file the indicates it in the data destination
    if simulated:
        with open(os.path.join(data_destination, "QASM"), 'w') as f: pass
    else:
        with open(os.path.join(data_destination, "ORIGINAL"), 'w') as f: pass

    # telling the user the the name of the directory that the data ended up in
    print("Done Copying data to:", data_destination)