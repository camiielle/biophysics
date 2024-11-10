import random

aa = "CMFILVWYAGTSQNEDHRKP"
hydrophobic_aa = aa[:10]
other_aa = aa[10:]

freq = [1.660, 2.370, 4.100, 5.810, 9.430, 6.580, 1.240, 3.190,
        7.580, 6.840, 5.670, 7.130, 3.970, 4.440, 6.360, 5.270,
        2.240, 5.160, 5.940, 4.920]
hydrophobic_freq=freq[:10]
other_freq=freq[10:]

def sequence_generator(length):
    sequence = []
    while len(sequence) < length:
        if len(sequence) % 2 == 0:
            sequence.append(random.choices(hydrophobic_aa, weights=hydrophobic_freq, k=1)[0])
        else:
            sequence.append(random.choices(other_aa, weights=other_freq, k=1)[0])
    sequence = ''.join(sequence)
    assert(len(sequence)==length)
    return sequence

length=500
print(sequence_generator(length))

