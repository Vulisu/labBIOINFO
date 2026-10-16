from collections import Counter
import re
import matplotlib.pyplot as plt

# Standard DNA to Amino Acid Genetic Code Table
DNA_TO_AA = {
    'TTT': 'F (Phenylalanine)', 'TTC': 'F (Phenylalanine)', 'TTA': 'L (Leucine)', 'TTG': 'L (Leucine)',
    'CTT': 'L (Leucine)', 'CTC': 'L (Leucine)', 'CTA': 'L (Leucine)', 'CTG': 'L (Leucine)',
    'ATT': 'I (Isoleucine)', 'ATC': 'I (Isoleucine)', 'ATA': 'I (Isoleucine)', 'ATG': 'M (Methionine)',
    'GTT': 'V (Valine)', 'GTC': 'V (Valine)', 'GTA': 'V (Valine)', 'GTG': 'V (Valine)',
    'TCT': 'S (Serine)', 'TCC': 'S (Serine)', 'TCA': 'S (Serine)', 'TCG': 'S (Serine)',
    'CCT': 'P (Proline)', 'CCC': 'P (Proline)', 'CCA': 'P (Proline)', 'CCG': 'P (Proline)',
    'ACT': 'T (Threonine)', 'ACC': 'T (Threonine)', 'ACA': 'T (Threonine)', 'ACG': 'T (Threonine)',
    'GCT': 'A (Alanine)', 'GCC': 'A (Alanine)', 'GCA': 'A (Alanine)', 'GCG': 'A (Alanine)',
    'TAT': 'Y (Tyrosine)', 'TAC': 'Y (Tyrosine)', 'TAA': '* (Stop)', 'TAG': '* (Stop)',
    'CAT': 'H (Histidine)', 'CAC': 'H (Histidine)', 'CAA': 'Q (Glutamine)', 'CAG': 'Q (Glutamine)',
    'AAT': 'N (Asparagine)', 'AAC': 'N (Asparagine)', 'AAA': 'K (Lysine)', 'AAG': 'K (Lysine)',
    'GAT': 'D (Aspartic Acid)', 'GAC': 'D (Aspartic Acid)', 'GAA': 'E (Glutamic Acid)', 'GAG': 'E (Glutamic Acid)',
    'TGT': 'C (Cysteine)', 'TGC': 'C (Cysteine)', 'TGA': '* (Stop)', 'TGG': 'W (Tryptophan)',
    'CGT': 'R (Arginine)', 'CGC': 'R (Arginine)', 'CGA': 'R (Arginine)', 'CGG': 'R (Arginine)',
    'AGT': 'S (Serine)', 'AGC': 'S (Serine)', 'AGA': 'R (Arginine)', 'AGG': 'R (Arginine)',
    'GGT': 'G (Glycine)', 'GGC': 'G (Glycine)', 'GGA': 'G (Glycine)', 'GGG': 'G (Glycine)'
}

def read_fasta(file_path):
    """Reads a FASTA/.fna file and returns the concatenated sequence in uppercase."""
    sequence_chunks = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence_chunks.append(line.upper())
    return "".join(sequence_chunks)

def analyze_codon_frequency(sequence):
    """Extracts codons using a non-overlapping step-3 sliding window."""
    sequence = re.sub(r'[^ATGC]', '', sequence)
    codons = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    return Counter(codons)

def main(fasta_file):
    sequence = read_fasta(fasta_file)
    codon_counts = analyze_codon_frequency(sequence)
    top_10 = codon_counts.most_common(10)
    
    # Prepare data for matplotlib
    codons = [item[0] for item in top_10]
    frequencies = [item[1] for item in top_10]
    
    # Create informative x-axis labels: DNA Codon / RNA Codon / Amino Acid Abbr
    x_labels = []
    for codon in codons:
        rna = codon.replace('T', 'U')
        aa_full = DNA_TO_AA.get(codon, "Unknown")
        aa_code = aa_full.split()[0] if aa_full != "Unknown" else "?"
        x_labels.append(f"{codon}\n({rna})\n[{aa_code}]")

    # Plotting configuration
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x_labels, frequencies, color='cornflowerblue', edgecolor='black', alpha=0.85)
    
    plt.xlabel('Codons (DNA / RNA / Amino Acid)', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency Count', fontsize=12, fontweight='bold')
    plt.title('Top 10 Most Frequent Codons in Genome', fontsize=14, fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Display the frequency value on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + (max(frequencies) * 0.01),
                 f'{height}',
                 ha='center', va='bottom', fontsize=10, fontweight='bold')
                 
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main("GCF_000855785.1_ViralMultiSegProj15044_genomic.fna")