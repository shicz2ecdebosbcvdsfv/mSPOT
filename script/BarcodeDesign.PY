import random
import csv

random.seed(42)

num_sequences = 112
seq_length = 24
min_levenshtein_distance = 10
max_homopolymer = 2
gc_content_range = (0.4, 0.6)

nucleotides = ['A', 'T', 'C', 'G']

def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)

    for i, c1 in enumerate(s1):
        current_row = [i + 1]

        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)

            current_row.append(
                min(insertions, deletions, substitutions)
            )

        previous_row = current_row

    return previous_row[-1]

def has_long_homopolymer(seq, max_run):
    count = 1
    last = seq[0]

    for nt in seq[1:]:
        if nt == last:
            count += 1

            if count > max_run:
                return True
        else:
            count = 1
            last = nt

    return False

def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)

def generate_random_sequence():
    return ''.join(
        random.choices(
            nucleotides,
            k=seq_length
        )
    )

sequences = []
attempts = 0
max_attempts = 100000

while len(sequences) < num_sequences and attempts < max_attempts:

    attempts += 1

    candidate = generate_random_sequence()

    if has_long_homopolymer(candidate, max_homopolymer):
        continue

    gc_frac = gc_content(candidate)

    if not (
        gc_content_range[0]
        <= gc_frac
        <= gc_content_range[1]
    ):
        continue

    if all(
        levenshtein(candidate, existing)
        >= min_levenshtein_distance
        for existing in sequences
    ):
        sequences.append(candidate)

output_file = 'dna_sequences_levenshtein.csv'

with open(output_file, 'w', newline='') as csvfile:

    writer = csv.writer(csvfile)

    writer.writerow(['Index', 'Sequence'])

    for idx, seq in enumerate(sequences, 1):
        writer.writerow([idx, seq])

print(
    f"Generation complete: "
    f"{len(sequences)} sequences saved to {output_file}"
)
