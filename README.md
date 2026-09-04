#mChIP-seq for Multiplex and Multifactorial Epigenomic Profiling Uncovers Cancer-specific Histone Features in Cellular and Circulating Nucleosomes

This repository contains the computational analyses associated with the study "mChIP-seq for Multiplex and Multifactorial Epigenomic Profiling Uncovers Cancer-specific Histone Features in Cellular and Circulating Nucleosomes"

mChIP-seq is a technology compatible with both cellular and cell-free samples for simultaneously profiling multifactorial epigenetic landscapes on multiple samples. To index samples, mChIP-seq uses the tailing and single-stranded ligation method to label the 3’ end of DNA in fragmented chromatin with sample-specific barcodes. The library structure includes barcode sequences on sample specificity, which are located at Read 2 for sample demultiplexing after PE sequencing.

[Uploading Library stucture.tif…]()

To demultiplex samples that are pooled for ChIP and library preparation, we can use existing software fastq-multx (https://github.com/brwnj/fastq-multx). After demultiplexing, each sample can be analyzed as traditional ChIP-seq data, such as bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for mapping, samtools (http://samtools.sourceforge.net/) for sam and bam file dealing, picard (http://broadinstitute.github.io/picard/) for PCR duplicates removing, macs2 (https://pypi.org/project/MACS2/) for peak calling.
