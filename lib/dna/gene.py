# genes.py - Wrappers around gene sequence functionality

# Standard library modules
import xml.etree.ElementTree as Tree

# Biopython modules
from Bio import Entrez
from Bio import SeqIO

# Library modules
from dna.cache import *
from dna.sequence import Sequence

class Gene:
    def defaults(key):
        return {
            'bacterial_strain': 'Ty2',
            'root': os.getcwd,
            'tmp': '/tmp/'
        }[key]

    def __init__(self, name, strain=defaults('bacterial_strain'), db='gene'):
        self.name = name
        self.description = name
        self.strain = strain
        self.db = db 
        self.cache = load_cache()
        
        self.gid = self.retrieve_gid()
        self.coordinates = self.retrieve_coordinates()
        self.sequence = self.retrieve_sequence()

        self.primers = {}
      
    def retrieve_gid(self):
        if self.name in cache['genes']:
            return cache['genes'][self.name]

        search = self.name + ' ' + self.strain
        handle = Entrez.esearch(db=self.db, term=search)
        record = Entrez.read(handle)
        
        # Our ID of interest should be the first result
        gid = record['IdList'] and record['IdList'][0]
        cache['genes'][self.name] = gid

        return gid

    def retrieve_coordinates(self):
        # We need an explicit gene ID or name, otherwise return
        if not self.gid:
            return

        elif self.gid in cache['coordinates']:
            return cache['coordinates'][self.gid]

        # XML identifiers for parsing Entrez file
        # 'Position' refers to location in the XML tree
        locus_position  = 'Entrezgene/Entrezgene_locus'
        region_position = 'Gene-commentary/Gene-commentary_seqs/' + \
                          'Seq-loc/Seq-loc_int/Seq-interval'
        start_position  = 'Seq-interval_from'
        end_position    = 'Seq-interval_to'
        gi_position     = 'Seq-interval_id/Seq-id/Seq-id_gi'
        strand_position = 'Seq-interval_strand/Na-strand'

        # Download and parse the correct region
        handle  = Entrez.efetch(db=self.db, id=self.gid, retmode='xml')
        locus   = Tree.parse(handle).getroot().find(locus_position)
        region  = locus.find(region_position)

        # Quick finder functions
        xml     = lambda p: region.find(p)
        val     = lambda p: xml(p).get('value')
        text    = lambda p: xml(p).text
        posn    = lambda p: int(text(p)) + 1

        # Find relevant information
        gi      = text(gi_position)
        strand  = val(strand_position)
        start   = posn(start_position)
        end     = posn(end_position)

        cache['coordinates'][self.gid] = (gi, strand, start, end)

        # Return a 4-tuple
        return gi, strand, start, end

    def retrieve_sequence(self):
        gi, strand, start, end = self.coordinates
        return str(Sequence(gi, strand, start, end))

    def retrieve_sequence_old(self):
        # Get the pertinent information from the parameters or
        # look them up based on the gene ID
        if not self.coordinates:
            return
              
        if self.coordinates in cache['sequences']:
            return cache['sequences'][self.coordinates]

        # Look up the sequence
        gi, strand, start, end = self.coordinates

        # Normalize the strands: plus = 1, minus=2
        strand = {'plus': 1, 'minus': 2}[strand.lower()]

        handle = Entrez.efetch(db='nucleotide', rettype='fasta', id=gi,
            seq_start=start, seq_stop=end, strand=strand)

        # Read, cache and return the sequence
        entry = SeqIO.read(handle, 'fasta')
        sequence = str(entry.seq)
        cache['sequences'][self.coordinates] = sequence

        return sequence

    def as_hash(self):
        gi, strand, start, end = self.coordinates
        gene = {
            'name': self.name,
            'gid': self.gid,
            'gi': gi,
            'strand': strand,
            'start': start,
            'end': end,
            'description': self.name,
            'sequence': self.sequence
        }

        for primer in self.primers.keys():
            gene[primer] = self.primers[primer]

        return gene
