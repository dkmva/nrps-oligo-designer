from cStringIO import StringIO
import os

from Bio import SeqIO, SearchIO
import argparse

import utils

# Constants

muscle_exe = 'lib/muscle'
infile1 = "lib/A_domains_muscle.fasta"
infile2 = "infile.fasta"
outfile = "lib/muscle.fasta"
reference_seq = "P0C062_A1"
startpos = 66
A_pos = [int(pos)-startpos for pos in "235 236 239 278 299 301 322 330 331".split()]

oligo_length = 90
min_homology = 15

hmmsearch = 'lib/hmmsearch'


class KnownCodes(object):

    def __init__(self, fasta):
        self.code_to_aa = {}
        self.aa_to_code = {}
        self.aa_to_dgn_code = {}

        self.parse_fasta(fasta)

    def parse_fasta(self, fasta):
        with open(fasta) as handle:
            for entry in SeqIO.parse(handle, 'fasta'):
                seq = str(entry.seq)
                aa = entry.id.split('__')[-1]

                self.code_to_aa[seq] = aa
                if aa not in self.aa_to_code:
                    self.aa_to_code[aa] = []
                    self.aa_to_dgn_code[aa] = []
                self.aa_to_code[aa].append(seq)
                self.aa_to_dgn_code[aa].append(utils.reverse_translate(seq))


class AdenylationDomain(object):

    def __init__(self, sequence, known=None, id='', revcom=False):
        self.nt_seq = sequence
        self.aa_seq = utils.translate_nucleotides(sequence)
        self.qseq = ''
        self.refseq = ''
        self.backwards = {}
        self.positions = []
        self.code = []
        self.specificity = 'unknown'
        self.library = {}
        self.oligos = {}
        self.id = id
        self.revcom = revcom

        self._muscle_align()
        self._get_positions()
        self._get_backwards()
        self._get_code()
        if known:
            self.get_specificity(known)

    def __repr__(self):
        return 'Adenylation domain {}: {}'.format(self.id, self.specificity)

    def _muscle_align(self):
        with open(infile2, 'w') as f:
            f.write('>qseq\n')
            f.write(self.aa_seq + '\n')

        os.system('{} -profile -quiet -in1 {} -in2 {} -out {}'.format(muscle_exe, infile1, infile2, outfile))

        with open(outfile) as handle:
            seqdict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            self.qseq = seqdict['qseq']
            self.refseq = seqdict[reference_seq]

        os.remove(infile2)
        return

    def _get_positions(self):
        a = 0
        b = 0
        for res in self.refseq:
            if a in A_pos and res !='-':
                self.positions.append(b)
            if res != '-':
                a += 1
            b += 1
        return

    def _get_backwards(self):
        c = 0
        for j in range(len(self.aa_seq)):
            while len(self.qseq) > j+c and self.qseq[j+c] == '-':
                c += 1
            self.backwards[j+c] = j
        return

    def _get_code(self):
        for pos in self.positions:
            res = self.aa_seq[self.backwards[pos]]
            start = self.backwards[pos]*3
            end = start + 3
            self.code.append((res, start, end))
        return

    def get_specificity(self, known):
        score = 0
        selected = ''
        stach = ''.join([r[0] for r in self.code]) + 'K'
        for code in known.code_to_aa:
            qscore = 0
            for q, r in  zip(code, stach):
                if q == r:
                    qscore += 1
            if qscore >= score:
                score = qscore
                selected = code

        self.specificity = known.code_to_aa[selected]
        return

    def make_library(self, known):

        if not self.specificity:
            self.get_specificity(known)

        for aa, dgns in known.aa_to_dgn_code.iteritems():
            if aa == self.specificity:
                continue
            codons = [self.nt_seq[res[1]:res[2]] for res in self.code]
            req_muations = 99
            req_codons = 99
            req_oligos = 99
            selected_dgn = ''
            cdn_chngs = []
            req_set = []

            for dgn in dgns:
                num_mutations = 0
                codon_changes = 0
                cdns = []
                for i, codon in enumerate(codons):
                    nt_changes = [nt in utils.dgn_to_nts[dg] for nt, dg in zip(codon, dgn[i])]
                    if False in nt_changes:
                        codon_changes += 1
                        cdns.append(True)
                    else:
                        cdns.append(False)
                    num_mutations += nt_changes.count(False)

                changes = [self.backwards[pos] for pos, tf in zip(self.positions, cdns) if tf][::-1]
                num_oligos = 1
                oligo_set = []
                fp = 0
                for i, pos in enumerate(changes):
                    if changes[fp]-pos > (oligo_length-2*min_homology)/3:
                        oligo_set.append(changes[fp:i][::-1])
                        fp = i
                        num_oligos += 1

                if not oligo_set or not changes[-1] in oligo_set[-1]:
                    oligo_set.append(changes[fp:][::-1])

                if num_oligos < req_oligos or (num_oligos == req_oligos and num_mutations < req_muations):
                    # TODO: Codon usage..
                    selected_dgn = dgn
                    req_muations = num_mutations
                    req_codons = codon_changes
                    cdn_chngs = cdns
                    req_oligos = num_oligos
                    req_set = oligo_set[::-1]

            self.library[aa] = {'dgn': selected_dgn,
                                'num_mutations': req_muations,
                                'codon_changes': req_codons,
                                'cdns': cdn_chngs,
                                'num_oligos': req_oligos,
                                'oligo_set': req_set}
        return

    def make_oligos(self):
        if not self.library:
            print 'Must have mutation library before making oligos.'
            return

        for aa in self.library:
            seq = self.nt_seq

            for s, e, c in [(res[1], res[2], dgn) for res, dgn, tf in zip(self.code, self.library[aa]['dgn'], self.library[aa]['cdns']) if tf]:
                # TODO: Codon usage..
                c =  ''.join([next(iter(utils.dgn_to_nts[nt])) for nt in c])
                seq = seq[:s] + c + seq[e:]

            oligo_positions = []
            mut_positions = []
            self.oligos[aa] = []
            for s in self.library[aa]['oligo_set']:
                mut_start = s[0]*3
                mut_end = s[-1]*3+3
                mut = seq[mut_start:mut_end]
                pre_len = (oligo_length - len(mut)) / 2

                oligo_start = mut_start-pre_len
                oligo_end = oligo_start + oligo_length

                oligo_seq = self.nt_seq[oligo_start:mut_start] + mut.lower() + self.nt_seq[mut_end:oligo_end]
                if self.revcom:
                    oligo_seq = utils.reverse_complement(oligo_seq)
                self.oligos[aa].append(oligo_seq)

                oligo_positions.append([oligo_start, oligo_end])
                mut_positions.append([mut_start, mut_end])

            for i, pos in enumerate(oligo_positions[1:]):
                if mut_positions[i][1] > pos[0]:
                    print 'Oligo clash detected..'

        return

    def write_oligo_string(self, string):
        for aa in self.oligos:
            for i, oligo in enumerate(self.oligos[aa], 1):
                string += '>{}:{}-{}:{}/{}\t{}\n'.format(self.id, self.specificity, aa, i, len(self.oligos[aa]), oligo)

        return string


def get_pepseq(dnaseq):
    pepseq = ''
    rf = 0
    for i in [-3,-2,-1,1,2,3]:
        tmpseq = utils.translate_nucleotides(dnaseq, reading_frame=i)
        if len(tmpseq) > len(pepseq):
            pepseq = tmpseq
            rf = i

    return pepseq, rf


def get_adenylation_domains(fasta, known=None, lagging_strand=False):
    adenylation_domains = []

    fasta_seqs = []
    for fs in SeqIO.parse(fasta, 'fasta'):
        revcom=False
        seq = str(fs.seq)
        pepseq, rf = get_pepseq(seq)
        if rf < 0 == lagging_strand:
            revcom=True
            seq = utils.reverse_complement(seq)
        fasta_seqs.append({'id': fs.id, 'seq': seq, 'pepseq': pepseq, 'rf': rf})
    for fs in fasta_seqs:
        utils.run_cmd([hmmsearch, '--domtblout', 'dump', os.path.abspath('lib/AMP-binding.hmm'), '-'],
                  '>header\n' + pepseq)
        with open('dump') as f:
            out = f.read()
        res_stream = StringIO(out)
        os.remove('dump')
        results = list(SearchIO.parse(res_stream, 'hmmsearch3-domtab'))

        for result in results:
            for i, hsp in enumerate(result.hsps, 1):
                s = hsp.hit_start
                e = hsp.hit_end

                adenylation_domains.append((AdenylationDomain(fs['seq'][s*3:e*3], known, '{}_{}'.format(fs['id'], i), revcom), s, e))

    return adenylation_domains

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('fasta', help='input sequence file')
    parser.add_argument('--lagging', action='store_true', help='lagging strand')

    args = parser.parse_args()

    known = 'lib/knowncodes.fasta'

    known = KnownCodes(known)

    domains = get_adenylation_domains(args.fasta, known, args.lagging)

    for dom in domains:
        dom[0].make_library(known)
        dom[0].make_oligos()
        print dom[0].write_oligo_string(''),
