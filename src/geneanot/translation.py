# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-locals
"""
Functions related to DNA and RNA translation.
"""
import random as rnd
import itertools as it
from collections import defaultdict
import numpy as np

# Basic DNA and RNA bp parameters
DNA_bases = {'T', 'C', 'A', 'G'}
RNA_bases = {s for s in DNA_bases if s != 'T'}.union('U')


def DNA2RNA(seq):
    """DNA to RNA """
    return seq.upper().replace('T', 'U')


def RNA2DNA(seq):
    """RNA to DNA """
    return seq.upper().replace('U', 'T')


def is_valid_DNA_seq(seq):
    """Returns True [False] is seq is a valid [not a valid] DNA sequence"""
    return set(seq.upper()) <= DNA_bases


def is_valid_RNA_seq(seq):
    """Returns True [False] is seq is a valid [not a valid] RNA sequence"""
    return set(seq.upper()) <= RNA_bases


def is_valid_bp_seq(seq, seq_type='DNA'):
    """Returns True [False] is seq is a valid [not a valid] seq_type sequence"""
    return is_valid_DNA_seq(seq) if seq_type=='DNA' else is_valid_RNA_seq(seq)


def reverse_complement(s: str, complement: dict | None = None) -> str:
    """Performs reverse-complement of a sequence. Default is a DNA sequence."""
    if complement is None:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    s_rev = s[::-1]
    lower = [b.islower() for b in list(s_rev)]
    bases = [complement.get(base, base) for base in list(s_rev.upper())]
    #rev_compl = ''.join([b.lower() if l else b for l, b in zip(lower, bases)])
    rev_compl = ''.join([bs.lower() if lw else bs for lw, bs in zip(lower, bases)])
    return rev_compl


def DNA_reverse_complement(s: str) -> str:
    """DNA reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return reverse_complement(s, complement=complement)


def RNA_reverse_complement(s: str) -> str:
    """RNA reverse complement"""
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return reverse_complement(s, complement=complement)


# Codon translation table
RNA_codon_table = {
    #    U             C             A             G
# U
    'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
    'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---', # UxA
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp', # UxG
# C
    'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
    'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
    'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
    'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG
# A
    'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
    'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG
# G
    'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
    'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
    'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}

# reverse translation ({aa : a codon of aa})
RNA_aa_table = {v: k for k, v in RNA_codon_table.items()}

# all 61 amin-acid-encoding codons
RNA_all_aa_codons = sorted([k for k, v in RNA_codon_table.items() if v != '---'])

# reverse translation with synonymous codons {aa : [all codons of aa]}
RNA_aa_synon_table = {v: [k for k, vk in RNA_codon_table.items() if vk == v] \
                      for v in set(RNA_codon_table.values())}

# AA code (Three-letter code : [IUPAC, Name])
AA_code = {
    'Ala': ['A', 'Alanine'],
    'Arg': ['R', 'Arginine'],
    'Asn': ['N', 'Asparagine'],
    'Asp': ['D', 'Aspartic-acid'],
    'Cys': ['C', 'Cysteine'],
    'Glu': ['E', 'Glutamic-acid'],
    'Gln': ['Q', 'Glutamine'],
    'Gly': ['G', 'Glycine'],
    'His': ['H', 'Histidine'],
    'Ile': ['I', 'Isoleucine'],
    'Leu': ['L', 'Leucine'],
    'Lys': ['K', 'Lysine'],
    'Met': ['M', 'Methionine'],
    'Phe': ['F', 'Phenylalanine'],
    'Pro': ['P', 'Proline'],
    'Ser': ['S', 'Serine'],
    'Thr': ['T', 'Threonine'],
    'Trp': ['W', 'Tryptophan'],
    'Tyr': ['Y', 'Tyrosine'],
    'Val': ['V', 'Valine']
}

# aa IUPAC to three-letter code (IUPAC AA: three letters AA)
AA_IUPAC23code = {v[0]: k for k, v in AA_code.items()}

# aa three-letter code to IUPAC (three letters AA: IUPAC AA)
AA_3code2IUPAC = {k: v[0] for k, v in AA_code.items()}

# IUPAC -> synonymous codon list
RNA_IUPAC_aa_synon_table = {AA_3code2IUPAC[k]: v for k, v in RNA_aa_synon_table.items() if k != '---'}


def convert_aa_IUPAC23(aa_seq):
    """Converts a IUPAC AA sequence to three-letter code AA sequence"""
    return ''.join([AA_IUPAC23code.get(i, '???') for i in list(aa_seq)])


def convert_aa_32IUPAC(aa_seq):
    """Converts a three-letter AA code sequence to IUPAC AA sequence"""
    return ''.join([AA_code.get(aa_seq[i:i+3], '?')[0] for i in range(0, len(aa_seq), 3)])


def translate_RNA_codon(codon):
    """Returns the amino acid corresponding to a codon"""
    return RNA_codon_table[DNA2RNA(codon)]


def translate_RNA_aa(aa):
    """Returns a codon corresponding to an animo acid"""
    return RNA_aa_table.get(aa, '???')


def translate(seq: str):
    """Returns the amino acid sequence corresponding to the RNA sequence seq"""
    translation = ''
    for n in range(0, len(seq) - (len(seq) % 3), 3): # every third base
        translation += translate_RNA_codon(seq[n:n+3])
    return translation

def translate_IUPAC(seq: str) -> str:
    """Returns the IUPAC amino acid sequence corresponding to the RNA sequence seq"""
    return convert_aa_32IUPAC(translate(seq))

# this can be used for both translation and reverse translation
def gen_translate(seq, func=translate_RNA_codon):
    """Performs translation based on a translation function"""
    translation = ''
    for n in range(0, len(seq) - (len(seq) % 3), 3): # every third base of three letters of AA
        translation += func(seq[n:n+3])
    return translation


def translate_in_frame(seq, framenum: int = 1):
    """Returns the translation of seq in a given reading frame"""
    #return translate(seq[framenum-1:])
    return gen_translate(seq[framenum-1:])


def print_translation_in_frame(seq, framenum: int = 1, prefix: str = ''):
    """Prints the translation of seq in reading frame framenum preceded by prefix"""
    print(prefix, framenum, ' ' * framenum, translate_in_frame(seq, framenum), sep='')


def print_translations(seq: str, prefix: str = '') -> None:
    """Prints the translation of seq in all three reading frames, each preceded by prefix"""
    print('\n',' ' * (len(prefix) + 2), seq, sep='', end=':\n')
    #[print_translation_in_frame(seq, fnum, prefix) for fnum in range(1,4)]
    for fnum in range(1,4):
        print_translation_in_frame(seq, fnum, prefix)

# Translation of ORF (start codon to stop codon)
def translate_with_open_reading_frames(seq: str, framenum: int = 1, open_enb: bool = False):
    """Returns the translation of seq in framenum (1, 2, or 3), with ---'s when not
    within an open reading frame"""
    translation = ''
    seqlength = len(seq) - (framenum - 1)
    for n in range(framenum-1, seqlength - (seqlength % 3), 3):
        codon = translate_RNA_codon(seq[n:n+3])
        open_enb = (open_enb or codon == "Met") and not (codon == "---")
        translation += codon if open_enb else "---"
    return translation


def print_translation_with_open_reading_frame(seq: str, framenum: int = 1, prefix: str = ''):
    """Prints ORF translation of a seq"""
    print(prefix, framenum, ' ' * framenum, translate_with_open_reading_frames(seq, framenum), sep='')


def print_translations_with_open_reading_frames(seq: str, prefix: str = ''):
    """Prints ORF translation of a seq for the three reading frame cases"""
    print('\n', ' ' * (len(prefix) + 2), seq, sep='')
    #[print_translation_with_open_reading_frame(seq, frame, prefix) for frame in range(1,4)]
    for frame in range(1,4):
        print_translation_with_open_reading_frame(seq, frame, prefix)
    
def gen_random_bp_seq(length: int, DNA: bool = True):
    """Generates a random sequence of DNA (default) or RNA nucleotides"""
    return ''.join(rnd.choices(list(DNA_bases), k=length)) if DNA else ''.join(rnd.choices(list(RNA_bases), k=length))


def gen_all_kmers(k: int, DNA: bool = True) -> list:
    """Generating all k-mers sequences given the value k."""
    bases = list(DNA_bases) if DNA else RNA_bases
    return ["".join(kmer) for kmer in list(it.product(bases, repeat=k))]


def gen_random_aa(length: int, IUPAC: bool = True):
    """Generates a random sequence of AA IUPAC letters (default) or a sequence of three letter-code AA"""
    return ''.join(rnd.choices(list(AA_IUPAC23code.keys()), k=length)) if IUPAC \
        else ''.join(rnd.choices(list(AA_code), k=length))


def seq_to_synon_codon_dist(aa_seq: list, codon_seq: list) -> tuple[dict, dict]:
    """
    Give a sequence of amino acids and the corresponding sequence of codons, the function
    computes the mapping between the unique amino-acids in the sequence and the codons used by the sequence.

    The function returns:
    1. aa_seq_code_prob: aa -> {codon: p}, where p is the probability of using codon to encode aa in the sequence.
    2. aa_seq_code: aa -> list of codons used to encode aa in the sequence (with repetitions).

    For example, assume that the amino-acid "I" appears 4 times in aa_seq and was
    encoded in codon_seq by 'AUU' (first time), 'AUU' (second time), 'ATA' (third time), and 'AUU' (forth time). Then:
    1. aa_seq_code['I'] = ['AUU', 'AUU', 'ATA', 'AUU']
    2. aa_seq_code_prob = {'AUU': 0.75, 'ATA': 0.25}

    aa_seq and codon_seq must not contain the stop codon.
    """
    # aa -> list of codon used (with repetitions)
    aa_seq_code = defaultdict(list)
    for a, c in zip(aa_seq, codon_seq):
        aa_seq_code[a].append(c)

    # converting to probabilities over all synonymous codons
    # aa -> {codon: prob}
    aa_seq_code_prob = {
        k: {cur_v: sum(np.array(v) == cur_v) / len(v) for cur_v in RNA_IUPAC_aa_synon_table[k]}
        for k, v in aa_seq_code.items()
    }
    return aa_seq_code_prob, dict(aa_seq_code)


def gen_rnd_codons_based_given_dict(aa_seq: list, aa_seq_code: dict) -> list:
    """
    Generates a random codon sequence (without the stop codon) that encodes
    aa_seq based on codon distribution given in aa_seq_code.

    aa_seq_code is the mapping table between an amino-acid and all allowed codons
    (with repetition) that encodes this amino acid. See the function seq_to_synon_codon_dist.

    aa_seq must not contain the stop codon.
    """
    seq = []
    for a in aa_seq:
        all_codons = aa_seq_code[a]  # list of codons (with repetition) that can be used to encode a
        seq.append(all_codons[np.random.randint(len(all_codons))])
    return seq


def seq_codon_distribution(codon_seq: list) -> tuple[dict, dict]:
    """
    Returns the distribution of all 61 codons based on the given codon sequence.

    codon_seq is a list of codons, where each codon is over the RNA alphabet (i.e. U is
    used and not T).
    """
    # make sure we work with RNA alphabet
    codon_rna_seq = [c.upper().replace('T', 'U') for c in codon_seq]

    # all 61 codon frequency
    codon_freq = {c: sum(np.array(codon_rna_seq) == c) for c in RNA_all_aa_codons}

    # normalized frequency
    codon_prob = {k: v / len(codon_rna_seq) for k, v in codon_freq.items()}

    return codon_freq, codon_prob
