# covid19_intrapatient_diversity
### Analysis of the population diversity of SARS-CoV-2 within and between individual patients

This data has the following potential applications:
* Genomic positions that do not mutate during infection may represent ideal targets for DNA-sequence-based therapeutics
* Protein residues that are under the highest levels of stabilizing selection could inform rational design of molecules to target these residues
* Comparing and contrasting how the virus evolves within a person versus globally within people gives insight into the viral lifestyle
* Estimating the number of viral particles acquired by a patient during infection has epidemiological implications

This analysis was performed primarily by Matt Olm (mattolm@stanford.edu) in Justin Sonnenburg's lab at Stanford University and Alex Crits-Christoph (crist-christoph@berkeley.edu) from Jill Banfield's lab at University of California, Berkeley.

**This work is on-going - please feel free to reach out with corrections, additions, or questions.**

# Data

Here we will link to all of the data that we generated so that other people and use / re-use as needed

### Parsed data-tables (available in the datatables folder)

* SRA_metadata_v1.csv - Metadata of the SRA samples used for intra-patient analysis

* COVID_genome_coverage_v1.csv - Genome-wide coverage, breadth, microsiversity, etc. of reads mapping to a COVID-19 reference genome

* COVID_positional_coverage_v1.csv - Position by position coverage and microdiversity along the genome

* COVID_SNPs_v1.csv - Locations of SNPs

# Methods

This method requires access to the raw sequencing reads generated when sequencing covid19 genomes. Usually these reads are used to generate a viral genome, the genome is deposited into a public database, and the raw reads are never uploaded publicly. **If you are involved in covid-19 genome sequencing efforts, please consider uploading the raw reads as well so that analyses like this can continue.**

### Sample identification

Here is where we'll put how I found the reads that I did

### Bioinformatics workflow

Here is where we'll put an overview of the methods used, the samples currently being analyzed

### Quality control

Some of the QC stuff that Matt has been doing.

# Intra- versus Inter- population diversity

This is really the only section that we'll fill out at first release. This will all be pi based

# Intra- versus Inter- population selection

This will essentially be repreating the above analysis with SNVs and dN/dS

# Estimating the number of viral particles transfered during transmission

Analysis on-going
