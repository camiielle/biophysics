# Write a program that goes through all consecutive (non-overlapping) triplets looking for
# stop codons. (Make sure you use the genetic code for DNA in the 5’-3’ direction.) Record
# the distance L between consecutive stop codons. Repeat this computation for the 3 different
# reading frames (0, +1, +2) in this direction. (You may skip calculations for the reverse
# strand, that is complementary to the given one and proceeding in the opposite direction.)

# Plot the distribution for the ORF lengths L calculated above, and compare it to that for random sequences

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
# NM: qui sto assumendo che la sequenza sia multiplo di 3 però
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

def find_distances(stop_positions):
    """Calculates the distances between consecutive stop codons in terms of codons."""
    distances = []
    for i in range(1, len(stop_positions)):
        distance = (stop_positions[i] - stop_positions[i - 1]-3) // 3  # Divide by 3 to get number of codons
        distances.append(distance)
    return distances

# generates a random sequence of nucleotides of length l
def generate_random_sequence(l):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(l))

file_path = './77177281caa9fedfeb5d8fdbfb33167d_U00096.fna'

real_sequence = read_file(file_path)
length = len(real_sequence)
real_distances_per_frame = []

# processes each of the three reading frames
# each frame produces a different set of codons
for frame in range(3):
    real_positions = find_stops(real_sequence, frame)   
    real_distances = find_distances(real_positions)
    real_distances_per_frame.append(real_distances) 
    
random_sequence = generate_random_sequence(length) 
random_positions = find_stops(random_sequence, 0)
random_distances = find_distances(random_positions)

def plot_scatter_from_probabilities(lengths, label, color):
    """Generates a scatter plot for ORF lengths and their probabilities."""
    length_counts = Counter(lengths)  # Count the frequency of each length
    total_orfs = sum(length_counts.values())  # Total number of ORFs for normalization
    length_probabilities = {length: count / total_orfs for length, count in length_counts.items()}  # Normalize to probabilities
    
    lengths_sorted = sorted(length_probabilities.items())  # Sort by ORF length
    lengths_list, probabilities = zip(*lengths_sorted)  # Unzip into two lists
    
    plt.scatter(lengths_list, probabilities, color=color, label=label, alpha=0.7, s=3)

# Calculate mean ORF lengths for each real genome frame
mean_frame_0 = statistics.mean(real_distances_per_frame[0])
mean_frame_1 = statistics.mean(real_distances_per_frame[1])
mean_frame_2 = statistics.mean(real_distances_per_frame[2])

# Calculate mean ORF length for random sequence (frame 0)
mean_random = statistics.mean(random_distances)

print(f"Mean ORF length for Frame 0: {mean_frame_0} codons")
print(f"Mean ORF length for Frame 1: {mean_frame_1} codons")
print(f"Mean ORF length for Frame 2: {mean_frame_2} codons")
print(f"Mean ORF length for Random Frame 0: {mean_random} codons")

colors = ['blue', 'green', 'orange']   
# Plotting
plt.figure(figsize=(10, 6))
    
# Plot ORF length vs frequency for each real frame
for frame in range(3):
    plot_scatter_from_probabilities(real_distances_per_frame[frame], label=f'Frame {frame}', color=colors[frame])
    
# Plot ORF length vs frequency for random frame 0
plot_scatter_from_probabilities(random_distances, label='Random', color='red')

# Add titles and labels
plt.title('ORFs length distribution', fontsize=14)
plt.xlabel('ORF length (codons)', fontsize=12)
plt.ylabel('Relative frequency', fontsize=12)
plt.legend()

# Set x-axis limit
# the mean is around 25 so i feel like this is a reasonable range
plt.xlim(-5, 180)
    
# Show the plot
plt.tight_layout()
plt.show()

# Function that calculates the fraction of ORFs greater than L_star
def fraction_orfs_greater_than_Lstar(lengths, L_star):
    total_orfs = len(lengths)
    if total_orfs == 0:
        return 0  # Avoid division by zero if there are no ORFs
    
    # Count the number of ORFs with length > L_star
    count_greater_than_Lstar = sum(1 for length in lengths if length > L_star)
    
    # Calculate the fraction
    fraction = count_greater_than_Lstar / total_orfs
    return fraction

# Set the range of L_star values
L_star_range = np.linspace(0, 200, 30) 
# Initialize lists to store fractions
random_fractions = []
real_fractions_per_frame = {0: [], 1: [], 2: []}  # Using a dictionary to store each frame

# Loop over L_star values and calculate fractions
for L_star in L_star_range:
    # Calculate fractions for random distances
    random_fractions.append(fraction_orfs_greater_than_Lstar(random_distances, L_star))
    
    # Calculate fractions for each real frame
    for i in range(3):
        real_fractions_per_frame[i].append(fraction_orfs_greater_than_Lstar(real_distances_per_frame[i], L_star))

# Plotting the results as scatter plots
plt.figure(figsize=(10, 6))

# Plot fractions for each frame as scatter plots
for i in range(3):
    plt.scatter(L_star_range, real_fractions_per_frame[i], label=f'Frame {i}', s=10)  # Adjust size (s) as needed

# Plot fractions for random data as a scatter plot
plt.scatter(L_star_range, random_fractions, label='Random', color='red', s=10)

# Set y-axis to log scale
plt.yscale('log')

# Set x-axis limit
plt.xlim(-5, 180)

# Add labels and title
plt.xlabel('L* (Codon Length)', fontsize=12)
plt.ylabel('Fraction of ORFs with L > L*', fontsize=12)
plt.title('Fraction of ORFs with L > L* as a Function of L* (Scatter)', fontsize=14)

# Add legend
plt.legend()

# Show plot
plt.grid(True)
plt.tight_layout()
plt.show()
