# Test file and backup file.
import random


def create_fake_genome(length):
    """ Creates a random sequence with a fake genome header"""
    fake_genome = ""
    with open("fake_genome.txt", "w") as writeto:
        writeto.write(">Fake_header_for_fake_genome" + "\n")
        for i in range(length):
            fake_genome += random.choice(["A", "C", "T", "G"])
            # Writes it to the new fake genome file in lines of length 50
        writeto.write('\n'.join(fake_genome[i:i + 50] for i in range(0, len(fake_genome), 50)))


create_fake_genome(30000)
