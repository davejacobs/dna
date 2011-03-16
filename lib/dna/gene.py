# gene.py - Wrapper around gene lookup functionality

# Standard library modules
import xml.etree.ElementTree as Tree

# Biopython modules
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

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
        self.upstream_sequence = {}
        self.downstream_sequence = {}

        self.primers = self.all_primers()
      
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

    def retrieve_upstream_sequence(self, length=50):
        # Look to see if we've already memoized this sequence
        # in our object
        if upstream = self.upstream_sequence[length]
            return upstream

        # Else look up the sequence per the normal lookup (which
        # will first try the cache and then ping Genbank
        gi, strand, start, end = self.coordinates
        upstream = self.upstream_sequence[length] = \
                str(Sequence(gi, strand, start-(length+1), start-1))

        return upstream

    def retrieve_downstream_sequence(self, length=50:
        # Look to see if we've already memoized this sequence
        # in our object
        if downstream = self.downstream_sequence[length]
            return downstream

        # Else look up the sequence per the normal lookup (which
        # will first try the cache and then ping Genbank
        gi, strand, start, end = self.coordinates
        downstream = self.downstream_sequence[length] = \
                str(Sequence(gi, strand, end+1, end+length+1))

        return downstream

    # Primer search and selection functions 
    # Constructs primers that will prime our plasmid (pKD4) and will have
    # extensions that are homologous with the regions of DNA that flank
    # gene.
    def homology_primer(self, relative_location='upstream', length=50, space=''):
        gi, strand, start, end = self.coordinates
        
        upstream = str(Sequence(gi, strand, start-length, start-1))
        downstream = str(Sequence(gi, strand, end+1, end+length))

        # Temporary solution to experimental parameters
        # (NOT final or how I would like it)
        experiment = {
            'forward_homology': 'GTGTAGGCTGGAGCTGCTTC',
            'reverse_homology': 'TAAGGAGGATATTCATATG'
        }

        # If we're on the negative strand, then our start and end
        # calculations were switched -- so switch the homologous
        # regions before adding any other sequences to them.
        if strand == 'minus':
            upstream, downstream = downstream, upstream

        upstream += space + experiment['forward_homology']
        downstream = str(Seq(experiment['reverse_homology'] + space + \
                downstream).reverse_complement())

        if relative_location == 'upstream':
            return upstream
        elif relative_location == 'downstream':
            return downstream
        else:
            return None

    # Returns the homology primers for a gene as a map
    def homology_primers(self, length=50):
        upstream = self.homology_primer('upstream', length)
        downstream = self.homology_primer('downstream', length)

        return {
            'homology-forward': upstream,
            'homology-reverse': downstream
        }

    # Returns internal primers for detecting the presence of gene
    # as a map
    def internal_primers(self):
        gi, strand, start, end = self.coordinates
        sequence = Sequence(gi, strand, start, end)

        return {
            'internal-forward': sequence.forward_primer(),
            'internal-reverse': sequence.reverse_primer()
        }

    # Returns an external primer for gene as a map. Distance
    # out from gene defaults to 500 bp
    def external_primers(self, distance=500):
        gi, strand, start, end = self.coordinates
        upstream = Sequence(gi, strand, start-distance, start-1)
        downstream = Sequence(gi, strand, end+1, end+distance)

        return {
            'external-forward': upstream.forward_primer(),
            'external-reverse': downstream.reverse_primer()
        }

    # Returns homology, internal and external primers for gene
    # as a map
    def all_primers(self):
        homology = self.homology_primers()
        internal = self.internal_primers()
        external = self.external_primers()

        # Why doesn't this work?
        # return reduce(dict.update, [homology, internal, external], {})

        primers = {}
        for p in [homology, internal, external]:
            primers.update(p)
        return primers

    # Returns hash representation of gene, so that we can interpolate
    # gene data into templates
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

        for primer in self.primers.items():
            gene[primer[0]] = primer[1]

        return gene
