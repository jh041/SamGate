# SamGate
Create input files for MegaSat

If you haven't noticed already, SamGate is an anagram of MegaSat, which is a software program designed to automate the calling of microsatellite alleles from Illumina fastq sequences (see Zhan et al. 2017. Molecular Ecology Resources 17:247-256).

SamGate takes as input the output files of MsatCommander--a classic program for detecting microsatellite loci from DNA sequence data (Faircloth 2008. Molecular Ecology Resources 8:92-94).

Before using SamGate, you will need to run MsatCommander on a fasta file reference (most likely a whole genome sequence) from which you will be harvesting microsatellite loci. MsatCommander is available here: https://github.com/brantfaircloth/msatcommander. Msatcommander will then generate two files: an msats file and a primer file.

To use SamGate, open a terminal window and run the python script:

$ python /Users/John/Documents/SamGate/SamGate.py

The script will prompt you for three files:

1) the original fasta file (.fa .fna .fasta)--the one that you used as input for Msat commander 
2) the msats file generated by msat commander
3) the primer file generated by msat commander

The result is a MegaSat formated file with the primers, motif and flanking regions for each locus:

Locus-Name	no_tail_left	Rev_Comp_of_no_tail_Right	ThreeFlank	FiveFlank	Repeat
46	ACCAACCTCAGTCATCAGGC	CATCTGCTTCCTGACTCCCA	ATTATTAAACTCAC	AGTATTGGTTCTGTTTTCATGTTTAGACGCTTCTCTTGATTTGTAA	ATC

Note: sometimes it is useful to concatenate MsatCommander output files before running SamGate. If you do this you will need to remove all
header lines before hand... in fact, I recommend removing all header and subheader lines before running SamGate, regardless. 
