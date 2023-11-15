from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC_skew

# Add your email address for NCBI's E-utilities
Entrez.email = 'rachelcasendra@gmail.com'

def fetch_genomic_data(accession_number):
    try:
        # Step 1: Fetch Genomic Data from NCBI
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Step 2: Reverse Sequences
        reversed_sequence = record.seq.reverse_complement()

        # Step 3: Calculate GC Skew
        gc_skew_values = GC_skew(reversed_sequence)

        return record, reversed_sequence, gc_skew_values

    except Exception as e:
        print(f"Error fetching genomic data: {e}")
        return None, None, None

def calculate_coverage(record):
    coverage = [0] * len(record.seq)

    for feature in record.features:
        # Assuming each feature contributes 1 to the coverage
        for position in feature.location:
            coverage[position] += 1

    return coverage

def plot_coverage(coverage, title="Coverage Plot"):
    plt.plot(range(len(coverage)), coverage)
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title(title)
    plt.show(block=True) 

def create_bed_file(record, gene_sequences, bed_file_path):
    with open(bed_file_path, 'w') as bed_file:
        for i, (gene_name, gene_sequence) in enumerate(gene_sequences):
            start = i
            end = i + len(gene_sequence)
            bed_file.write(f"{record.id}\t{start}\t{end}\t{gene_name}\n")

def create_track_file(track_file_path):
    with open(track_file_path, 'w') as track_file:
        track_file.write("track exampleTrack\n")
        track_file.write("type bed 3\n")
        track_file.write("shortLabel Your Custom Track\n")
        track_file.write("longLabel Your Custom Track Description\n")
        track_file.write("visibility hide\n")       

def visualize_gc_skew(record, reversed_sequence, gc_skew_values):
    # Step 4: Visualize GC Skew
    plt.plot(range(len(gc_skew_values)), gc_skew_values)
    plt.xlabel("Position")
    plt.ylabel("GC Skew")
    plt.title(f"GC Skew of {record.id}")

    # Save the plot to a file
    plt.savefig("gc_skew_plot.png")

def extract_gene_sequences(record):
    gene_sequences = []
    
    for feature in record.features:
        if feature.type == 'gene':
            gene_sequence = feature.extract(record.seq)
            gene_sequences.append((feature.qualifiers['gene'][0], gene_sequence))
    
    return gene_sequences

def calculate_gene_lengths(gene_sequences):
    gene_lengths = {gene_name: len(sequence) for gene_name, sequence in gene_sequences}
    return gene_lengths

def visualize_gene_locations(record, gene_sequences):
    fig, ax = plt.subplots()
    
    for i, (gene_name, gene_sequence) in enumerate(gene_sequences):
        ax.text(i, len(record.seq) + 10, gene_name, rotation=45, ha='right', va='bottom')
        ax.add_patch(plt.Rectangle((i, 0), len(gene_sequence), len(record.seq), color='lightgray'))

    ax.set_xlim(0, len(record.seq))
    ax.set_ylim(0, len(record.seq) + 20)
    ax.set_xlabel("Position")
    ax.set_title(f"Gene Locations in {record.id}")

    plt.savefig("gene_locations_plot.png")   

def annotate_genes(record):
    # Step 5: Annotate Genes
    gene_features = [feature for feature in record.features if feature.type == 'gene']

    for gene_feature in gene_features:
        print(f"Gene: {gene_feature.qualifiers['gene'][0]}")
        print(f"Location: {gene_feature.location}")
        print(f"Description: {gene_feature.qualifiers.get('product',['No description available'])[0]}\n")

if __name__ == "__main__":
    # Replace 'your_accession_number' with the actual accession number you want to analyze
    accession_number = "JX573431.1"

    # Fetch genomic data, reverse sequences, and calculate GC skew
    record, reversed_sequence, gc_skew_values = fetch_genomic_data(accession_number)

    if record:
        # Visualize GC skew and save the plot to a file
        visualize_gc_skew(record, reversed_sequence, gc_skew_values)

        # Annotate and print information about genes
        annotate_genes(record)

        # Extract gene sequences
        gene_sequences = extract_gene_sequences(record)

        # Create a BED file for gene locations
        create_bed_file(record, gene_sequences, "gene_locations.bed")

        # Create a track file for the UCSC Genome Browser
        create_track_file("trackDb.txt")

        # Calculate coverage
        coverage = calculate_coverage(record)

        # Plot coverage
        plot_coverage(coverage, title=f"Coverage Plot of {record.id}")
