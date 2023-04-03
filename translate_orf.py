import sys
import re

def vet_nucleotide_sequence(sequence):
    rna_pattern_str = r'^[AaUuCcGg]*$'
    dna_pattern_str = r'^[AaTtCcGg]*$'

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))

def vet_codon(codon):
    codon_pattern_str = r'^[AaUuCcGg]{3}$'

    codon_pattern = re.compile(codon_pattern_str)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))

def translate_codon(codon):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    return genetic_code[codon]

def find_first_orf(sequence, start_codons=['AUG'], stop_codons=['UAA', 'UAG', 'UGA']):
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    # Get copies of everything in uppercase
    upper_sequence = sequence.upper()
    upper_start_codons = [codon.upper() for codon in start_codons]
    upper_stop_codons = [codon.upper() for codon in stop_codons]

    orfs = []

    #Find all start codons
    for start_codon in upper_start_codons:
        start_indexes = [m.start() for m in re.finditer(start_codon, upper_sequence)]

    # Find all stop codons after each start codon
    for start_index in start_indexes:
        found_stop_codon = False
        for i in range(start_index + 3, len(upper_sequence), 3):
            codon = upper_sequence[i:i+3]
            if codon in upper_stop_codons:
                found_stop_codon = True
                orfs.append(upper_sequence[start_index:i+3])
                break

        # If no stop codon was found, assume the ORF goes to the end of the sequence
        if not found_stop_codon:
            orfs.append(upper_sequence[start_index:])

    #Return the longest ORF
    if orfs:
        return max(orfs, key=len)
    else:
        return None

def main():
    sequence = input("Enter a nucleotide sequence: ")
    orf = find_first_orf(sequence)
    if orf:
        protein_sequence = ''.join([translate_codon(orf[i:i+3]) for i in range(0, len(orf), 3)])
        print("Protein sequence:", protein_sequence)
    else:
        print("No ORF found.")

if __name__ == '__main__':
    main()

