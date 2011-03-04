from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import SeqRecord
from Bio.Emboss.Applications import Primer3Commandline
from dna.cache import *
from dna.defaults import *

class Sequence:
    def __init__(self, gi, strand, start, end):
        self.gi = gi
        self.strand = strand
        self.start = start
        self.end = end
        self.coordinates = (gi, strand, start, end)

        self.sequence = self.retrieve_sequence()

    def __str__(self):
        return self.sequence

    def retrieve_sequence(self):
        # Get the pertinent information from the parameters or
        # look them up based on the gene ID
        if self.coordinates in cache['sequences']:
            return cache['sequences'][self.coordinates]

        # Normalize the strands: plus = 1, minus=2
        self.strand = {'plus': 1, 'minus': 2}[self.strand.lower()]

        handle = Entrez.efetch(db='nucleotide', rettype='fasta', id=self.gi,
            seq_start=self.start, seq_stop=self.end, strand=self.strand)

        # Read, cache and return the sequence
        entry = SeqIO.read(handle, 'fasta')
        sequence = str(entry.seq)
        cache['sequences'][self.coordinates] = sequence

        return sequence

    def primer(self, direction='forward'):
        # Stub method because eprimer3 is not available on my local machine
        # at the moment

        return 'stub-sequence'

        tmp_file = defaults['tmp'] + self.gi + '-' + \
                self.strand + str(self.start) + str(self.end)

        primer3_input_file = tmp_file
        primer3_output_file = tmp_file + '-output'

        record = SeqRecord(Seq(self.sequence), id='', description='')

        with open(primer3_input_file, 'w') as f:
            SeqIO.write(record, f, 'fasta')

        primercl = Primer3Commandline(sequence=primer3_input_file, auto=True,
                hybridprobe=True)

        primercl.osizeopt = 20
        primercl.psizeopt = 200
        primercl.outfile = primer3_output_file

        stdout, stderr = primercl()

        unpack = lambda r: [[p.forward_seq, p.reverse_seq] for p in r.primers]
        
        with open(primer3_output_file, 'r') as f:
            record = Primer3.parse(f).next()
            all_primers = [[p.forward_seq, p.reverse_seq] for p in record.primers]

        if direction == 'forward':
            return all_primers[0][0]

        elif direction == 'reverse':
            return all_primers[0][1]

        else:
            return None

    def forward_primer(self):
        return self.primer('forward')

    def reverse_primer(self):
        return self.primer('reverse')
