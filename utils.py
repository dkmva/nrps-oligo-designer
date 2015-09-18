import subprocess
from string import maketrans

translation_table = {'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', 'ctt': 'L',
                     'ctc': 'L', 'cta': 'L', 'ctg': 'L', 'att': 'I', 'atc': 'I',
                     'ata': 'I', 'atg': 'M', 'gtt': 'V', 'gtc': 'V', 'gta': 'V',
                     'gtg': 'V', 'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
                     'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 'act': 'T',
                     'acc': 'T', 'aca': 'T', 'acg': 'T', 'gct': 'A', 'gcc': 'A',
                     'gca': 'A', 'gcg': 'A', 'tat': 'Y', 'tac': 'Y', 'taa': 'X',
                     'tag': 'X', 'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q',
                     'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K', 'gat': 'D',
                     'gac': 'D', 'gaa': 'E', 'gag': 'E', 'tgt': 'C', 'tgc': 'C',
                     'tga': 'X', 'tgg': 'W', 'cgt': 'R', 'cgc': 'R', 'cga': 'R',
                     'cgg': 'R', 'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R',
                     'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'}

degenerate_nucleotides = (
    ("A", frozenset(["A"])),
    ("T", frozenset(["T"])),
    ("G", frozenset(["G"])),
    ("C", frozenset(["C"])),
    ("R", frozenset(["A", "G"])),
    ("M", frozenset(["A", "C"])),
    ("Y", frozenset(["C", "T"])),
    ("K", frozenset(["G", "T"])),
    ("S", frozenset(["C", "G"])),
    ("W", frozenset(["A", "T"])),
    ("B", frozenset(["C", "G", "T"])),
    ("D", frozenset(["A", "G", "T"])),
    ("H", frozenset(["A", "C", "T"])),
    ("V", frozenset(["A", "C", "G"])),
    ("N", frozenset(["A", "C", "G", "T"]))
)

nts_to_dgn = {nts: dgn for dgn, nts in degenerate_nucleotides}
dgn_to_nts = {dgn: nts for dgn, nts in degenerate_nucleotides}

default_dgn_table = {'A': 'GCN', 'C': 'TGY', 'E': 'GAR', '$': 'TRR',
                     'G': 'GGN', 'F': 'TTY', 'I': 'ATH', 'H': 'CAY',
                     'K': 'AAR', 'M': 'ATG', 'L': 'YTN', 'N': 'AAY',
                     'Q': 'CAR', 'P': 'CCN', 'S': 'WSN', 'R': 'MGN',
                     'T': 'ACN', 'W': 'TGG', 'V': 'GTN', 'Y': 'TAY',
                     'D': 'GAY'}


def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
        raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
    return out


def complement(x):
    transtable = maketrans('ATGCRYMKBVDHatgcrymkbvdh',
                           'TACGYRKMVBHDtacgyrkmvbhd')
    return str(x).translate(transtable)


def reverse(x):
    return str(x)[::-1]


def reverse_complement(x):

    x = complement(x)
    x = reverse(x)

    return x


def translate_nucleotides(sequence, reading_frame=1, stop=True):
    sequence = str(sequence).lower()
    peptide_sequence = ''

    if reading_frame < 0:
        sequence = reverse_complement(sequence)
        reading_frame = abs(reading_frame)

    for i in range((reading_frame - 1), len(sequence), 3):
        codon = sequence[i:(i+3)]
        if len(codon) < 3:
            break
        try:
            residue = translation_table[codon]
            if stop and residue == 'X':
                break
            peptide_sequence += residue
        except KeyError:
            raise KeyError('Unknown codon: "{}"'.format(codon))

    return peptide_sequence


def reverse_translate(seq):
    return [default_dgn_table[AA] for AA in seq]