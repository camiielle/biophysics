# EXERCISE 2:
# a) Write a program that goes through all consecutive (non-overlapping) triplets looking for
# stop codons. (Make sure you use the genetic code for DNA in the 5’-3’ direction.) Record
# the distance L between consecutive stop codons. Repeat this computation for the 3 different
# reading frames (0, +1, +2) in this direction. (You may skip calculations for the reverse
# strand, that is complementary to the given one and proceeding in the opposite direction.)

# b) Plot the distribution for the ORF lengths L calculated above, and compare it to that for random sequences

import random
import matplotlib.pyplot as plt
from collections import Counter
import statistics
import numpy as np

# takes a .fna file and reads the sequence of the E.coli genome


def read_file(file_path):
    sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

# finds the stop codons when given a specified frame
# I take into account the circularity


def find_stops(sequence, frame):
    stop_codons = ["TAA", "TGA", "TAG"]
    stop_positions = []
    for i in range(frame, len(sequence) + frame - 2, 3):
        if i + 3 <= len(sequence):
            codon = sequence[i:i+3]
        else:
            codon = sequence[i:] + sequence[:frame]
        if codon in stop_codons:
            stop_positions.append(i)

    return stop_positions

# calculates distances between consecutive stop codons (in terms of codons, not bases)


def find_distances(stop_positions):
    distances = []
    for i in range(1, len(stop_positions)):
        # I divide by 3 to get distance in terms of codons
        distance = (stop_positions[i] - stop_positions[i - 1] - 3) // 3
        distances.append(distance)
    return distances


def generate_random_sequence(l):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(l))


def plot_ORF_length_distribution(lengths, label, color):
    length_counts = Counter(lengths)
    N_ORFs = sum(length_counts.values())
    length_frequencies = {length: count /
                          N_ORFs for length, count in length_counts.items()}
    lengths_sorted = sorted(length_frequencies.items())
    lengths_list, frequencies = zip(*lengths_sorted)
    plt.scatter(lengths_list, frequencies, color=color, label=label, s=3)


file_path = './77177281caa9fedfeb5d8fdbfb33167d_U00096.fna'
real_sequence = read_file(file_path)
length = len(real_sequence)

# processes each of the three reading frames
# each frame produces a different set of codons
real_distances_per_frame = []
for frame in range(3):
    real_positions = find_stops(real_sequence, frame)
    real_distances = find_distances(real_positions)
    real_distances_per_frame.append(real_distances)

random_sequence = generate_random_sequence(length)
random_positions = find_stops(random_sequence, 0)
random_distances = find_distances(random_positions)

# calculates mean ORF length
mean_frame_0 = statistics.mean(real_distances_per_frame[0])
mean_frame_1 = statistics.mean(real_distances_per_frame[1])
mean_frame_2 = statistics.mean(real_distances_per_frame[2])
mean_random = statistics.mean(random_distances)

# plotting
plt.figure(figsize=(10, 6))
colors = ['blue', 'green', 'orange']
for frame in range(3):
    plot_ORF_length_distribution(
        real_distances_per_frame[frame], label=f'Frame {frame}', color=colors[frame])

plot_ORF_length_distribution(random_distances, label='Random', color='red')

mean_text = f"Mean ORF length (codons):\nFrame 0: {mean_frame_0:.2f}\nFrame 1: {mean_frame_1:.2f}\nFrame 2: {mean_frame_2:.2f}\nRandom: {mean_random:.2f}"
plt.text(0.77, 0.15, mean_text, transform=plt.gca().transAxes, fontsize=10,
         bbox=dict(facecolor='white', alpha=0.5))

plt.title('ORFs Length Distribution', fontsize=16)
plt.xlabel('ORF length (codons)', fontsize=11)
plt.ylabel('Relative frequency', fontsize=11)
plt.legend()
# the mean is around 25 so i feel like this is a reasonable range
plt.xlim(-5, 180)
plt.grid(True)
plt.tight_layout()
plt.show()

# c) Estimate a cut-off value Lcut, above which the ORFs are statistically significant, i.e. the number of observed
# ORFs with L > Lcut is much greater than expected by chance.

# calculates the fraction of ORFs greater than L_cut


def fraction_greater_than_Lcut(lengths, L_cut):
    N_ORFs = len(lengths)
    if N_ORFs == 0:
        return 0  # i want to avoid division by zero
    count = 0
    for length in lengths:
        if length > L_cut:
            count += 1
    fraction = count / N_ORFs
    return fraction


L_cut_range = np.linspace(0, 200, 30)
random_fractions = []
real_fractions_per_frame = {0: [], 1: [], 2: []}

for L_cut in L_cut_range:
    random_fractions.append(
        fraction_greater_than_Lcut(random_distances, L_cut))
    for i in range(3):
        real_fractions_per_frame[i].append(
            fraction_greater_than_Lcut(real_distances_per_frame[i], L_cut))

# plotting
plt.figure(figsize=(10, 6))
for i in range(3):
    plt.scatter(
        L_cut_range, real_fractions_per_frame[i], label=f'Frame {i}', s=10)
plt.scatter(L_cut_range, random_fractions, label='Random', color='red', s=10)
# log scale so it's easier to graphically see the deviation from random case
plt.yscale('log')
plt.xlim(-5, 180)  # same range as plot before
plt.xlabel('L_cut (codons)', fontsize=11)
plt.ylabel('Relative fraction of ORFs with L > L_cut', fontsize=11)
plt.title('ORFs Statistical Significance', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

# plot for calculating L_cut with alternative method using function found in exercise 1
p_s = 3/64


def random_function(L):
    return (1-p_s)**L * p_s


plt.figure(figsize=(10, 6))
for frame in range(3):
    plot_ORF_length_distribution(
        real_distances_per_frame[frame], label=f'Frame {frame}', color=colors[frame])
L_range = np.arange(0, 180)
random_values = [random_function(L) for L in L_range]
plt.plot(L_range, random_values, color='red', label='Random Function')
plt.title('ORFs Statistical Significance (alternative method)', fontsize=16)
plt.xlabel('ORF length (codons)', fontsize=11)
plt.ylabel('Relative frequency', fontsize=11)
plt.yscale('log')
plt.legend()
plt.xlim(-5, 180)  # same range as plot before
plt.grid(True)
plt.minorticks_on()
plt.grid(True, which='minor', linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.show()

# EXERCISE 4:
# please note that sum functions, variables may repeat and be re-defined (e.g. read_file, sequence)
# because i wrote the code and used it as 2 separate files for the two exercises.
# If you want to try and run it i recommend separating them.

# a) Calculate the frequencies of the four nucleotides in the genome.

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
    spaced_pair_list.extend([sequence[-i] + sequence[n-i]
                            for i in range(1, n+1)])  # circularity
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
    matrix = construct_correlation_matrix(
        frequencies, spaced_pair_frequencies(sequence, n))
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
