# SARS-CoV-2 viral genomics

The purpose of this GitHub page is to share data, analyses, and results related to the analysis of genomic variation in the SARS-CoV-2 genome. **This work is on-going - please feel free to reach out with corrections, additions, or questions.**

This analysis was performed primarily by Matt Olm (<mattolm@stanford.edu>) in Justin Sonnenburg's lab at Stanford University and Alex Crits-Christoph (<crits-christoph@berkeley.edu>) in Jill Banfield's lab at University of California, Berkeley and analysis is broken up into the following major sections:

## [Summary](#Summary-1)

An overview of our results and links to quickly download some of the data parsed and generated during this study, including a sequence alignment of publicly available genomes on NCBI, locations of genomic SNPs, a table of gene annotations, and a table SRA run locations of raw reads.

## [Results](#Results-1)

Analysis based on comparing the covid19 genomes assembled from different patients (interpatient variation) and based on comparing the viral genetic diversity found within single individuals during infection (intrapatient variation).

## [Data availability](#Data-availability-1)

Links to the data generated during this study for others to use

## [Methods](#Methods-1)

Detailed methods on how analysis was performed.

## [Citations and acknowledgements](#citations-and-acknowledgements-1)

This work completely depends on the scientists and universities that originally sequenced these genomes and made their data publicly available.

Note: Intrapatient analyses require access to the raw sequencing reads generated when sequencing covid19 genomes. Usually these reads are used to generate a consensus viral genome, the genome is deposited into a public database, and the raw reads are never uploaded publicly. **If you are involved in covid-19 genome sequencing efforts, please consider uploading the raw reads as well so that analyses like this can continue.**

# Summary

This will be filled in last

# Results

## Introduction

As viruses replicate within their hosts during infection, they quickly mutate into genetically diverse populations. This is especially true for viruses with RNA genomes (as is the case with SARS-CoV-2). These populations are usually summarized as a single genome representing the consensus genome from a particular patient, and analysis of **inter-patient variation** (differences among the consensus genomes from different patients) is useful for understanding how the outbreak spreads, and how the virus evolves globally.

The variation among individuals in single host population is referred to as **intra-patient variation**, intraspecfic variation, or microdiversity. These type of variation is less commonly studied, but analysis of this data has the following potential applications:

* Identification of the genomic loci least likely to mutate during infection. Could be useful for designing universal primers / probes.

* Comparison of viral evolution within individuals versus global evolution. This can be useful for understanding how the viral evolutionary pressures and function.

* Estimation of the number of viral particles acquired to start infection and quantifying genetic diversity transferred during transmission.

## Inter-patient SARS-CoV-2 genome variation

[interFigure1]: results/interpatient_pi.png

To analyze inter-patient variation, we downloaded publicly available SARS-CoV-2 genomes, processed and filtered them, and generated a multiple sequence alignment. Details can be found in the [Methods](#Methods-1) second, and the genomes and alignment can be downloaded from the [Data availability](#Data-availability-1) section.

After genome alignment, the following cites and genes have the highest [nucleotide diversity (π)](https://en.wikipedia.org/wiki/Nucleotide_diversity) in the alignment:

![interFigure1][interFigure1]

The figure on the right shows that some positions along the genome are very conserved (low nucleotide diversity), and that other sites are much more variable. Sites could be variable because there's a lower selective cost to mutating them (genetic drift) or because there's a selective advantage to mutating them.

The first on the left shows that nucleotide diversity varies greatly between the different proteins produced by the genome. Proteins with the name "nsp" are "non-structural proteins" that make up the RNA-depended RNA polymerase and the rest of the replisome. The "step_loop" regions and the protein "nsp11" have much higher nucleotide diversity than other regions, and the spike protein "S" (a common vaccine target) has relatively low nucleotide diversity.  

## Intra-patient SARS-CoV-2 nucleotide diversity

[Figure1]: results/QC_boxplots_v2.png
[Figure2.1]: results/CoverageDistrubution_Center_ID_v2.png
[Figure2.2]: results/MicrodiversityDistrubution_Center_ID_v2.png
[Figure3]: results/CodonMicrodiversity_min_50_v2.png
[Figure3.2]: results/CodonMicrodiversity_min_1000_v2.png
[Figure4]: results/GeneMicrodiversity_RankOrder_v2.png
[Figure4.2]: results/GeneMicrodiversity_RankOrder_Coverage_v2.png
[Figure5]: results/GeneMicrodiversity_RankOrder_interIntra_v2.png
[Figure5.2]: results/GeneMicrodiversity_scatter_interIntra_v2.png

### Effect of sequencing protocol on resulting reads

The reads from the SRA were prepared using a number of different RNA extraction methods, library preparation methods, and DNA sequencing machines. To identify potential biases associated with different methods, we plotted the distribution of mapping quality metrics as compared to library preparation metadata (as retrieved from the SRA).

![Figure 1][Figure1]

The figure shows that the sequencing protocol does have an impact on metrics such as nucleotide diversity and coverage standard deviation, but it is difficult to disentangle the effects of the different variables. Overall it does not look like any particular method has a systemic bias associated with it.

We also plotted out the nucleotide [coverage distribution][Figure2.1] and [nucleotide diversity distribution][Figure2.2]. This data led us to conclude that ~1000x coverage is needed for a smooth nucleotide diversity distribution, and that ~50x coverage is needed for a smooth coverage distribution. Going forward we restricted our analysis primarily to the 24 samples that have >= 50x coverage overall.

### Intra-patient nucleotide diversity

To investigate which genomic positions have the highest intra-patient diversity (microdiversity), we use a previously established [measure of nucleotide diversity (π)](https://en.wikipedia.org/wiki/Nucleotide_diversity) (Nei and Li 1979). This metric is calculated for each position along the genome using inStrain, has a several advantages over typical SNP-based analyses, including avoiding SNP-calling thresholds that differ from program to program and allowing consideration of every read/position in the read mapping.

As a check to ensure everything is working correctly, we compared the nucleotide diversity of each codon position:

![Figure 3][Figure3]

The height of each bar represents the average normalized nucleotide diversity of all bases with that codon position among all genes and all samples with >=50x coverage, and black bars show the 95% confidence intervals. All codon positions have significantly different nucleotide diversities from one another (Wilcoxon rank-sum test; see notebook for details).

As expected, nucleotide diversity is highest at third codon position. This is because [mutations to the third codon position are much less likely to change the resulting amino acid than changes to the first or second position](https://genomevolution.org/wiki/index.php/Codon_wobble_positions), and thus, mutations to the third codon position less likely to have deleterious effects on the virus. The fact that our pipeline is able to detect this signal shows that the evolutionary pressures faced by the virus during infection are reflected in the resulting sequencing microdiversity. [A similar figure was also made for samples at over 1000x coverage][Figure3.2] that gave similar results.

We next looked at the overall microdiversity of each gene in the genome:

![Figure 4][Figure4]

For each set of reads analyzed we ranked each annotation by its intra-patient nucleotide diversity. The annotation with the lowest nucleotide diversity was ranked 1, second lowest 2, ect., all the way up to 31 (the total number of annotations). We then plotted the average rank of each annotation. Black bars represent 95% confidence intervals.  

[A second plot in which coverage is plotted along with nucleotide diversity][Figure4.2] is available as well to asses it as a potential confounding factor.

## Comparison of intra- and inter- patient diversity

![Figure 5][Figure5]

![Figure 5.2][Figure5.2]

## Conclusions

A couple of bullet points







# Data availability

[NCBI_genomes]: data/interpatient/ncbi_mar20.fna
[genome_alignment]: data/interpatient/mar20_filtered.aln


[SRA_metadata]: data/datatables/intrapatient/SRA_metadata_v1.csv
[COVID_genome_coverage]: data/datatables/intrapatient/COVID_genome_coverage_v2.csv
[COVID_positional_coverage]: data/datatables/intrapatient/COVID_positional_coverage_v2.csv.gz

[COVID_genes_table]: data/reference_genome/COVID_genes_table_v2.tsv
[COVID_genes_positional]: data/reference_genome/COVID_genes_positional_v2.csv

[interpatient_nucleotide_diversity]: data/interpatient/nucleotide_diversity.txt
[interpatient_snps]: data/interpatient/interpatient_snps.txt



## Raw data

* [NCBI_genomes][NCBI_genomes] - Fasta file of genomes used in this study

* [genome_alignment][genome_alignment] - Alignment of COVID-19 genomes used in this study

## Datatables


### Reference genome annotation

* [COVID_genes_table][COVID_genes_table] - Locations of gene-annotations on reference genome. Parsed from the reference GenBank file. **0-based indexing**

* [COVID_genes_positional][COVID_genes_positional] - Each position in the genome and all associated annotations. Parsed from the [COVID_genes_table][COVID_genes_table]. Providing because parsing is made difficult due to the confusing gene structure

### Inter-patient data

* [SRA_metadata][SRA_metadata] - Metadata of the SRA samples used for intra-patient analysis

* [COVID_genome_coverage][COVID_genome_coverage] - Genome-wide coverage, breadth, nucleotide diversity, etc. of reads mapping to a COVID-19 reference genome

* [COVID_positional_coverage][COVID_positional_coverage] - Position by position coverage and nucleotide diversity along the genome in sequencing reads

### Inter-patient data

* [interpatient_nucleotide_diversity][interpatient_nucleotide_diversity] - Nucleotide diversity of every position in the genome alignment

* [interpatient_snps][interpatient_snps] - Describes the position, nucleotides, and frequencies of each substitution across NCBI genomes.

## Data analysis notebooks

# Methods

## Inter-patient variation

### Downloading genomes

We are working with viral genomes from NCBI's data hub:
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049

But a similar analysis could be done on the larger collection of sequences available in [GISAID](gisaid.org). Currently there are about ~170 genomes in NCBI and over 800 in GISAID.

A FASTA file of March 20th currently available NCBI sequences can be found in [/data/interpatient/ncbi_mar20.fna][NCBI_genomes]

### Filtering

The script `./data/interpatient/filter_seqs.py` will perform some quality checks on the NCBI sequences. Firstly, several of the sequences are too short, and do not represent full viral genomes. We filter out all sequences < 29 Kb, as the viral genome is ~30 Kb. There are also several sequences that are too short in GISAID.

This script then filters out viral genomes that are too divergent to represent accurate genomes from the 2019 SARS-CoV-2 pandemic. This second step requires the program `fastANI` be installed in your path. It compares all genomes to the reference SARS-CoV-2 (`./data/interpatient/reference.fna`) and removes those that are <99% nucleotide identity.

### Creating a genome alignment

We then align the resulting filtered sequences using MAFFT - `mafft --thread 16 mar20_filtered_seqs.fna > mar20_filtered.aln` to generate a full genome alignment of the high quality viral genomes.

We generate a phylogenetic tree of these sequences using IQTree:

`iqtree -nt 10 -s mar20_filtered.aln`
and the generated tree is available in `./data/interpatient/mar20_filtered.iqtree`.

**The viral genome alignment is made available in [./data/interpatient/mar20_filtered.aln][genome_alignment]**

### Calculating metrics from the genome alignment

We conduct a number of analyses on the alignment. In general, we want to know which positions in the genome are variable, and how variable they are. The notebook `./data/interpatient/interpatient_pi.ipynb` calculates nucleotide diversity in the alignment, with respect to the reference genome's gff annotation. It produces the figures illustrating nucleotide diversity across the entire genome and by each gene and/or coding sequence.

Tabular data on the raw nucleotide diversity per reference position can be found in `./data/interpatient/nucleotide_diversity.txt`. Note that this calculation ignores any insertions or deletions from the reference and only tracks single nucleotide polymorphisms.

The notebook `./data/interpatient/interpatient_snps.ipynb` then generated the following table of genomic substitutions across patients in `./`:

```
ref_pos	nucleotide	frequency
1	A	0.9333333333333333
1	C	0.05333333333333334
1	T	0.013333333333333334
2	C	0.01282051282051282
2	T	0.9871794871794872
3	C	0.012658227848101266
3	T	0.9746835443037974
```
Describing the position, nucleotides, and frequencies of each substitution across NCBI genomes.

## Intra-patient variation

A typical pipeline for SARS-CoV-2 genome sequencing involves first generating a large number of DNA sequencing reads (short sequences of DNA that come from a single viral particle) and then assembling this data into a patient consensus genome (a representation of the most common viral genotype in a sample). Consensus genomes can be used for identifying outbreak clusters and global spread of the virus, but to understand the microdiversity present in a single sample, we use the raw DNA sequencing reads. Microdiversity (intrapatient variation) was profiled using the program [inStrain](https://github.com/MrOlm/instrain)

### Downloading and processing raw reads

Raw reads were identified by searching the term "SARS-CoV-2" in the NCBI SRA, and selecting the platform "Illumina". This was last performed on **3/17/2020** and led to the identification of **36** runs. This can be accessed at the following [URL](https://www.ncbi.nlm.nih.gov/sra/?term=(%22Severe+acute+respiratory+syndrome+coronavirus+2%22%5BOrganism%5D+OR+SARS-CoV-2%5BAll+Fields%5D)+AND+%22platform+illumina%22%5BProperties%5D). This information was downloaded into a .csv file by clicking "Send to:" -> Choose Destination: "File" -> Format "RunInfo" -> "Create File".

SRA files were downloaded using the URLs in the above file using the command:

```
$ cat SraRunInfo.csv | tail -n +2 | awk -F ',' '{print "wget " $10}' | bash
```

SRA files were converted into .fastq files using the command:

```
$ for s in $(ls SRR* | grep -v SRR11140750); do fastq-dump $s --split-files; done
```

For paired reads (reads where files ending in `_1.fastq` and `_2.fastq` were both created), the BBtools command `repair.sh` was used to remove un-paired reads using the command:

```
repair.sh in={fastq1} in2={fastq2} out={fastq1.repair} out2={fastq2.repair}
```

All reads were then trimmed and filtered using the BBtools command `bbduk.sh`:

```
# For paired reads:
repair.sh in={fastq1.repair} in2={fastq2.repair} out={fastq1.bbduk} out2={fastq2.bbduk}

# For unpaired reads:
bbduk.sh in={fastq1} out={fastq1.bbduk}
```

A table listing all SRA files currently analyzed in this study is available at [datatables/SRA_metadata_v1.csv](datatables/SRA_metadata_v1.csv)

### Mapping reads to a reference genome and calculating microdiversity

Microdiversity is calculated using the program [inStrain](https://github.com/MrOlm/instrain) which uses a `.bam` file as its input. A `.bam` file describes where each individual read or read-pair best maps to reference genome.

The reference genome used in this study is [NCBI Reference Sequence: NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). The `.fasta` file was downloaded for mapping by clicking "Sent to:" -> "Complete Record" -> Choose Destination: "File" -> Format "FASTA" -> "Create File". A `Genbank` format version listing gene locations was also downloaded using the same process but clicking Format "GenBank (full)".

Reads were mapped to the reference genome using the program Bowtie2 with the following commands:

```
# Make an index for mapping
$ bowtie2-build NC_045512.2.fasta NC_045512.2.fasta.bt2

# Example command for paired reads
$ bowtie2 -x /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.fasta.bt2 -1 /home/mattolm/user_data/Covid_19/reads/filtered/SRR11278092_bbduk.fastq --no-unal -S /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11278092.sam -p 10 2> /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11278092.sam.log

# Example command for unpaired reads
$ bowtie2 -x /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.fasta.bt2 -U /home/mattolm/user_data/Covid_19/reads/filtered/SRR11278092_bbduk.fastq --no-unal -S /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11278092.sam -p 10 2> /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11278092.sam.log
```

InStrain was then run on the resulting `.sam` files using the following commands (inStrain automatically converts `.sam` for `.bam` files):

```
# Example for paired reads
inStrain profile /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11241255.sam /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.fasta -o /home/mattolm/user_data/Covid_19/inStrain/profiles/NC_045512.2.fasta.bt2-vs-SRR11241255.sam.IS -p 10 --pairing_filter non_discordant -g /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.gb --skip_mm_profiling

# Example for unpaired reads
inStrain profile /home/mattolm/user_data/Covid_19/inStrain/mapping/NC_045512.2.fasta.bt2-vs-SRR11278167.sam /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.fasta -o /home/mattolm/user_data/Covid_19/inStrain/profiles/NC_045512.2.fasta.bt2-vs-SRR11278167.sam.IS -p 10 --pairing_filter non_discordant -g /home/mattolm/user_data/Covid_19/genomes/NC_045512.2.gb --skip_mm_profiling
```

### Processing inStrain results

Here's how I did that with links to jupyter notebooks. Include all that weird genes stuff.








# Citations and acknowledgements

This work completely depends on the scientists that originally sequenced these genomes and made their data publicly available. For the intra-patient analysis, this includes:

* Technological University of Pereira

* University of Washington

* The Wuhan Institute of Virology

* Beijing Institute of Genomics

* The University of Wisconsin–Madison
