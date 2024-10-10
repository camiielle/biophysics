# EXERCISE 4:
# a) Calculate the frequencies of the four nucleotides in the genome.

import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

# takes a .fna file and reads the sequence of the E.coli genome
def read_file(file_path):
    sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

# calculates nucleotides' relative frequencies
def nucleotide_frequencies(sequence):
    nucleotide_counter = Counter(sequence)
    N_nucleotides = sum(nucleotide_counter.values())
    return {nucleotide: count / N_nucleotides for nucleotide, count in nucleotide_counter.items()}

file_path = './77177281caa9fedfeb5d8fdbfb33167d_U00096.fna'
sequence = read_file(file_path)
frequencies = nucleotide_frequencies(sequence)
for nucleotide, freq in frequencies.items():
    print(f"Nucleotide: {nucleotide}, Frequency: {freq:.4f}")


# b) Write a program to count all 16 possible pairs of neighboring bases (e.g AT); hence obtain (1) the joint probabilities pXY
# and construct the 4×4 matrix of correlations cXY = pXY/(pXpY).

# calculates nucleotides' relative frequencies
def pair_frequencies(sequence):
    pair_list = [sequence[i:i+2] for i in range(len(sequence) - 1)]
    pair_list.append(sequence[-1] + sequence[0])  # circularity
    pair_counter = Counter(pair_list)
    assert len(pair_counter) == 16  # i expect 16 possible pairs
    N_pairs = len(sequence)
    return {pair: count / N_pairs for pair, count in pair_counter.items()}


def construct_correlation_matrix(nucleotide_freq, pair_freq):
    nucleotides = ['A', 'G', 'C', 'T']
    correlation_matrix = np.zeros((4, 4))
    for i, n1 in enumerate(nucleotides):
        for j, n2 in enumerate(nucleotides):
            pair = n1 + n2
            pXY = pair_freq.get(pair)
            pX = nucleotide_freq.get(n1)
            pY = nucleotide_freq.get(n2)
            correlation_matrix[i, j] = pXY / (pX * pY)
    return correlation_matrix


def construct_symbolic_matrix():
    nucleotides = ['A', 'G', 'C', 'T']
    matrix = np.full((4, 4), '', dtype='<U2')
    for i, n1 in enumerate(nucleotides):
        for j, n2 in enumerate(nucleotides):
            matrix[i, j] = n1+n2
    return matrix

# prints the matrix
print("\n4x4 Correlation Matrix:\n")
print(construct_symbolic_matrix())
print()
np.set_printoptions(precision=4)
print(construct_correlation_matrix(
    frequencies, pair_frequencies(sequence)))


# (c) Repeat the above calculation for nucleotides that are further neighbors, and find the corresponding matrices
# (e.g. consider next nearest neighbor locations j and j + 2 to calculate cXY(2)).
# How do correlations decay as a function of the separation n?

# calculates n-spaced pair frequencies
def spaced_pair_frequencies(sequence, n):
    spaced_pair_list = [sequence[i] + sequence[i+n]
                        for i in range(len(sequence) - n)]
    spaced_pair_list.extend([sequence[-i] + sequence[n-i] for i in range(1,n+1)]) # circularity
    spaced_pair_counts = Counter(spaced_pair_list)
    assert len(spaced_pair_counts) == 16
    N_spaced_pairs = len(sequence)
    return {pair: count / N_spaced_pairs for pair, count in spaced_pair_counts.items()}

# plotting 
y_axis_AA = []
y_axis_AG = []
y_axis_AC = []
y_axis_AT = []
for n in range(1, 41):
    matrix = construct_correlation_matrix(frequencies, spaced_pair_frequencies(sequence, n))
    y_axis_AA.append(matrix[0, 0])  
    y_axis_AG.append(matrix[0, 1]) 
    y_axis_AC.append(matrix[0, 2]) 
    y_axis_AT.append(matrix[0, 3]) 

x_axis = list(range(1, 41))
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

axs[0, 0].plot(x_axis, y_axis_AA, '-o', color='blue', label='AA')
axs[0, 0].set_title('AA')
axs[0, 0].grid(True)
axs[0, 0].set_ylabel('Correlation')  

axs[0, 1].plot(x_axis, y_axis_AG, '-o', color='orange', label='AG')
axs[0, 1].set_title('AG')
axs[0, 1].grid(True)

axs[1, 0].plot(x_axis, y_axis_AC, '-o', color='red', label='AC')
axs[1, 0].set_title('AC')
axs[1, 0].grid(True)
axs[1, 0].set_xlabel('Separation distance n') 
axs[1, 0].set_ylabel('Correlation')  

axs[1, 1].plot(x_axis, y_axis_AT, '-o', color='green', label='AT')
axs[1, 1].set_title('AT')
axs[1, 1].grid(True)
axs[1, 1].set_xlabel('Separation distance n')  

fig.suptitle("Correlation as a Function of Separation Distance", fontsize=16)
plt.tight_layout(rect=[0, 0, 0.95, 1])  
plt.show()
