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
    for i in range (frame, len(sequence) + frame - 2, 3):
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
        distance = (stop_positions[i] - stop_positions[i - 1]- 3) // 3  # I divide by 3 to get distance in terms of codons
        distances.append(distance)
    return distances

def generate_random_sequence(l):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(l))

def plot_ORF_length_distribution(lengths, label, color):
    length_counts = Counter(lengths)  # counts the frequency of each length
    N_ORFs = sum(length_counts.values())  
    length_frequencies = {length: count / N_ORFs for length, count in length_counts.items()}  
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

# Calculate mean ORF length 
mean_frame_0 = statistics.mean(real_distances_per_frame[0])
mean_frame_1 = statistics.mean(real_distances_per_frame[1])
mean_frame_2 = statistics.mean(real_distances_per_frame[2])
mean_random = statistics.mean(random_distances)

# Creating the plot
plt.figure(figsize=(10, 6))
colors = ['blue', 'green', 'orange']  
for frame in range(3):
    plot_ORF_length_distribution(real_distances_per_frame[frame], label=f'Frame {frame}', color=colors[frame])

plot_ORF_length_distribution(random_distances, label='Random', color='red')

mean_text = f"Mean ORF length (codons):\nFrame 0: {mean_frame_0:.2f}\nFrame 1: {mean_frame_1:.2f}\nFrame 2: {mean_frame_2:.2f}\nRandom: {mean_random:.2f}"
plt.text(0.77, 0.15, mean_text, transform=plt.gca().transAxes, fontsize=10,
         bbox=dict(facecolor='white', alpha=0.5))

plt.title('ORFs Length Distribution', fontsize=16)
plt.xlabel('ORF length (codons)', fontsize=11)
plt.ylabel('Relative frequency', fontsize=11)
plt.legend()
plt.xlim(-5, 180) # the mean is around 25 so i feel like this is a reasonable range
plt.grid(True)
plt.tight_layout()
plt.show()



# c) Estimate a cut-off value Lcut, above which the ORFs are statistically significant, i.e. the number of observed 
# ORFs with L > Lcut is much greater than expected by chance.

# function that calculates the fraction of ORFs greater than L_cut
def fraction_greater_than_Lcut(lengths, L_cut):
    N_ORFs = len(lengths)
    if N_ORFs == 0:
        return 0  # i want to avoid division by zero 
    count = 0
    for length in lengths:
        if length > L_cut:
            count+=1
    fraction = count / N_ORFs
    return fraction

L_cut_range = np.linspace(0, 200, 30) 
random_fractions = []
real_fractions_per_frame = {0: [], 1: [], 2: []}  

for L_cut in L_cut_range:
    random_fractions.append(fraction_greater_than_Lcut(random_distances, L_cut))
    for i in range(3):
        real_fractions_per_frame[i].append(fraction_greater_than_Lcut(real_distances_per_frame[i], L_cut))

# plotting 
plt.figure(figsize=(10, 6))
for i in range(3):
    plt.scatter(L_cut_range, real_fractions_per_frame[i], label=f'Frame {i}', s=10)  
plt.scatter(L_cut_range, random_fractions, label='Random', color='red', s=10)
plt.yscale('log') # log scale so it's easier to graphically see the deviation from random case
plt.xlim(-5, 180) # same range as plot before
plt.xlabel('L_cut (codons)', fontsize=11)
plt.ylabel('Relative fraction of ORFs with L > L_cut', fontsize=11)
plt.title('ORFs Statistical Significance', fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
