
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


file_path = './77177281caa9fedfeb5d8fdbfb33167d_U00096.fna'
sequence = read_file(file_path)


def count_codons_starting_with_AA_AG_AT_AC(sequence):
    count_AA = 0
    count_AG = 0
    count_AT = 0
    count_AC = 0

    for i in range(0, len(sequence) - 2, 3):  
        codon = sequence[i:i+3]

        if codon.startswith('AA'):
            count_AA += 1
        elif codon.startswith('AG'):
            count_AG += 1
        elif codon.startswith('AT'):
            count_AT += 1
        elif codon.startswith('AC'):
            count_AC += 1

    total = count_AA+count_AG+count_AC+count_AT
    return {
        'AA': count_AA/total,
        'AG': count_AG/total,
        'AT': count_AT/total,
        'AC': count_AC/total
    }

result = count_codons_starting_with_AA_AG_AT_AC(sequence)

print(result)

#reading frame one
sequence_one = sequence[1:] + sequence[0]
sequence_2 = sequence[2:] + sequence[0:2]

result_one = count_codons_starting_with_AA_AG_AT_AC(sequence_one)
result_2 = count_codons_starting_with_AA_AG_AT_AC(sequence_2)

print(result_one)
print(result_2)