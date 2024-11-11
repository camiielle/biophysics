import random
import numpy as np

# code for exercise 5b)

aa = "CMFILVWYAGTSQNEDHRKP"
hydrophobic_aa = aa[:10]
polar_aa = aa[10:]

freq = [1.660, 2.370, 4.100, 5.810, 9.430, 6.580, 1.240, 3.190,
        7.580, 6.840, 5.670, 7.130, 3.970, 4.440, 6.360, 5.270,
        2.240, 5.160, 5.940, 4.920]
hydrophobic_freq=freq[:10]
polar_freq=freq[10:]

def sequence_generator(length):
    sequence = []
    while len(sequence) < length:
        if len(sequence) % 2 == 0:
            sequence.append(random.choices(hydrophobic_aa, weights=hydrophobic_freq, k=1)[0])
        else:
            sequence.append(random.choices(polar_aa, weights=polar_freq, k=1)[0])
    sequence = ''.join(sequence)
    assert(len(sequence)==length)
    return sequence

length=1500
print(sequence_generator(length))


# code for exercise 5c)

h2b = "MTDKITKKKRNETYSIYIYKVLRQVHPKIGVSSKAMNIMNSFVNDLFERLVSESYNLSNSSRSKTLTAREIQTSVRLVIPGELAKHSVSEGTKAVAKYRSSI"
total_aa = len(h2b)
positions = range(total_aa)

def percentage_mutation(X, h2b):
    number_of_mutations = int(total_aa * X / 100)
    mutation_positions = random.sample(positions, k=number_of_mutations)
    h2b_list = list(h2b)
    
    for i in mutation_positions:
        h2b_list[i] = random.choice(aa)  
    
    return ''.join(h2b_list)

x_values = np.linspace(0, 100, 11, dtype=int)
for x in x_values:
    mutated_sequence = percentage_mutation(x, h2b)
    print(f"{x}% mutations: {mutated_sequence}")
