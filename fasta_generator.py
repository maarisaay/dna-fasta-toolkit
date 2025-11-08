import random
import json
from datetime import datetime

def generate_random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def calculate_stats(seq):
    length = len(seq)
    counts = {nuc: seq.count(nuc) for nuc in 'ACGT'}
    percentages = {nuc: (counts[nuc] / length) * 100 for nuc in 'ACGT'}
    cg = counts['C'] + counts['G']
    at = counts['A'] + counts['T']
    cg_at_ratio = cg / at if at != 0 else 0
    return counts, percentages, cg_at_ratio

# formatowanie FASTA
def wrap_fasta(seq, width):
    return ''.join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    batch_count = int(input("Ile sekwencji wygenerować? "))

    fmt = input("Wybierz format zapisu (f = FASTA, j = JSON): ").strip().lower()
    if fmt == 'j':
        output_file = "multi_sequences.json"
    else:
        output_file = "multi_sequences.fasta"

    metadata = f"##Wygenerowano: {datetime.now().isoformat()} | Tool: Python"

    entries = []

    if fmt != 'j':
        fasta = open(output_file, 'w')
        fasta.write(metadata + "")

    all_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_length = 0

    for i in range(1, batch_count + 1):
        print(f"Sekwencja {i}/{batch_count}")

        length = int(input('Podaj długość sekwencji: '))
        seq_id = input('Podaj ID sekwencji: ')
        description = input('Podaj opis sekwencji: ')
        name = input('Podaj imię: ')

        while True:
            dna_seq = generate_random_dna(length)

            # oblicz GC
            counts_tmp, _, _ = calculate_stats(dna_seq)
            gc_tmp = (counts_tmp['G'] + counts_tmp['C']) / length * 100

            # sprawdź powtórzenia > 4
            long_repeat = any(n*5 in dna_seq for n in "ACGT")

            warnings = []
            if gc_tmp < 30 or gc_tmp > 70:
                warnings.append(f"GC% = {gc_tmp:.2f}% (poza zakresem 30–70%)")
            if long_repeat:
                warnings.append("Wykryto powtórzenia nukleotydu dłuższe niż 4.")

            if warnings:
                print("\n OSTRZEŻENIA Quality Control:")
                for w in warnings:
                    print(" -", w)
                retry = input("Wygenerować sekwencję ponownie? (t/n): ").strip().lower()
                if retry == 't':
                    continue
            break
        insert_position = random.randint(0, length)
        dna_with_name = dna_seq[:insert_position] + name + dna_seq[insert_position:]

        counts, percentages, cg_at_ratio = calculate_stats(dna_seq)

        wrapped = wrap_fasta(dna_with_name, 60)

        if fmt != 'j':
            fasta.write(f">{seq_id} {description}")
            fasta.write(wrapped + "")

        entries.append({
            'id': seq_id,
            'description': description,
            'sequence': dna_with_name,
            'wrapped': wrapped,
            'insert_position': insert_position,
            'stats': {
                'counts': counts,
                'percentages': percentages,
                'cg_at': cg_at_ratio
            }
        })

        for n in 'ACGT':
            all_counts[n] += counts[n]
        total_length += length

        print(f"Pozycja wstawienia imienia: {insert_position}")
        for nuc, pct in percentages.items():
            print(f"{nuc}: {pct:.1f}%")
        print(f"CG/AT: {cg_at_ratio:.2f}")

    print("Statystyki")
    for n in 'ACGT':
        print(f"{n}: {all_counts[n] / total_length * 100:.2f}%")
    cg = all_counts['C'] + all_counts['G']
    at = all_counts['A'] + all_counts['T']
    print(f"CG/AT: {cg / at if at != 0 else 0:.2f}")

    if fmt == 'j':
        with open(output_file, 'w') as jf:
            json.dump({
                'metadata': metadata,
                'entries': entries
            }, jf, indent=4)
    else:
        fasta.close()

    print(f"Wszystkie sekwencje zapisano do pliku: {output_file}")

if __name__ == "__main__":
    main()
