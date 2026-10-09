CODON_TABLE = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
    'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',
}

def translate_rna(sequence: str) -> str:

    clean_seq = sequence.upper().replace('T', 'U').replace(' ', '').replace('\n', '')
    
    start_index = clean_seq.find('AUG')
    if start_index == -1:
        return "Error: No start codon (AUG) found."

    protein = ""

    for i in range(start_index, len(clean_seq) - 2, 3):
        codon = clean_seq[i:i+3]
        amino_acid = CODON_TABLE.get(codon, '?')

   
        if amino_acid == '*':
            break

        protein += amino_acid

    return protein



user_sequence = input("Enter your RNA or DNA sequence: ")
result = translate_rna(user_sequence)

print("\nResulting Amino Acid Sequence:")
print(result)