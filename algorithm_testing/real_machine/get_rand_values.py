import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import random
import os
import re

def plot(title, omegas, results):
    # figure that displays information on the omega database
    plt.figure(title, tight_layout=True, figsize=(12, 5))

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


with open(os.path.join("..", "..", "generated_data", "database")) as f:
    nums = re.split('[ \n]', f.read()) # splits the file's contents at newlines and spaces

nums = list(filter(("").__ne__, nums)) # cleans up the data
for i, item in enumerate(nums):
    nums[i] = item.split("->")
    nums[i] = list(map(float, nums[i]))

omegas = np.array(nums)[:, 0]
results = np.array(nums)[:, 1]

input_omegas = []
input_results = []

for i in range(1000):
    index = random.randint(0, len(omegas))
    input_omegas.append(omegas[index])
    input_results.append(results[index])

with open("omegas", "w") as f:
    for omega in input_omegas:
        f.write(str(omega) + "\n")

plot("Random Omega", input_omegas, input_results)

plot("Source Omega", omegas, results)

plt.show()