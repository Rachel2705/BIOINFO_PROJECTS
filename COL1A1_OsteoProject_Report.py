import Bio
print(Bio.__version__)

from io import StringIO
from Bio import Entrez

# Provide your email address to NCBI
Entrez.email = "rachelcasendra@gmail.com"

# Fetch the sequence using the gene ID
gene_id = "1277"
handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
record = handle.read()
handle.close()

print(record)

from Bio import SeqIO

# Parse the GenBank record
genbank_record = SeqIO.read(StringIO(record), "genbank")

from Bio.Seq import Seq
from Bio.SeqUtils import seq3

# Access information like sequence, features, etc.
sequence = genbank_record.seq
features = genbank_record.features

print("Sequence:", sequence)
print("Features:", features)

from Bio.SeqUtils import GC, translate

# Translate DNA to protein
protein_sequence = Seq(str(sequence)).translate()
print("Protein Sequence:", seq3(protein_sequence))



# Calculate GC content of the original DNA sequence
gc_content_dna = Bio.SeqUtils.GC_fraction(sequence)
print("GC Content of DNA Sequence:", gc_content_dna)


from Bio.Blast import NCBIWWW, NCBIXML

# Perform a BLAST search
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
blast_record = NCBIXML.read(result_handle)

# Access BLAST results
for alignment in blast_record.alignments:
    print("Alignment Title:", alignment.title)

from Bio.PDB import PDBList, PDBParser

# Download a PDB file
pdb_id = "1a3n"
pdb_list = PDBList()
pdb_file = pdb_list.retrieve_pdb_file(pdb_id)

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure(pdb_id, pdb_file)

# Access information about the structure
for model in structure:
    for chain in model:
        for residue in chain:
            print("Residue:", residue)

# Structure Analysis (Optional):
from Bio.PDB import PDBList, PDBParser

# Download a PDB file
pdb_id = "1a3n"
pdb_list = PDBList()
pdb_file = pdb_list.retrieve_pdb_file(pdb_id)

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure(pdb_id, pdb_file)

# Access information about the structure
for model in structure:
    for chain in model:
        for residue in chain:
            print("Residue:", residue)
         
