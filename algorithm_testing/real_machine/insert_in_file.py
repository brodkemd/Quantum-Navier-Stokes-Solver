import re

with open("input_omegas", "r") as f: nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces
nums = list(filter(("").__ne__, nums)) # cleans up the data

# gets the lines of script file
with open("algo.py", 'r') as f: lines = f.readlines()

to_insert = ", ".join(nums)
lines[0] = f"nums = [{to_insert}]\n"

# writes the lines to the script file
with open("algo.py", 'w') as f: f.writelines(lines)
