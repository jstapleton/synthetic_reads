# This makefile will reproduce the Gelsemium sythetic read assembly.
# Typically, the pairedDict.txt and contig_list.txt steps are run separately as queued jobs on a computer cluster.
# qsub scripts for these two steps are provided as HPCC_Hasher.sub and HPCC_Spades.sub
# Prior to running HPCC_Spades.sub, run:
#	mkdir splitDicts
#	cd splitDicts
#	python ~/source/dictSplitter.py ../unpairedDict.txt 900
#	cd ..



all: forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq unpairedDict.txt contig_list.txt contig_lengths


clean:
	-rm -f barcodePairs.txt confirmedPairs.txt forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq unpairedDict.txt contig_list.txt contig_lengths

forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Jasmine_E_CATTTT_L002_R1_001.fastq.gz Jasmine_E_CATTTT_L002_R2_001.fastq.gz forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired1.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
	rm reverse_unpaired.fq

forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Jasmine-E_ATCACG_L002_R1_001.fastq.gz Jasmine-E_ATCACG_L002_R2_001.fastq.gz forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired2.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

forward_paired_merged.fq: forward_paired1.fq forward_paired2.fq
	cat forward_paired1.fq forward_paired2.fq > forward_paired_merged.fq

reverse_paired_merged.fq: reverse_paired1.fq reverse_paired2.fq
	cat reverse_paired1.fq reverse_paired2.fq > reverse_paired_merged.fq

forward_unpaired_merged.fq: forward_unpaired1.fq forward_unpaired2.fq
	cat forward_unpaired1.fq forward_unpaired2.fq > forward_unpaired_merged.fq

unpairedDict.txt: forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq
	python barcodeHasher.py forward_paired_merged.fq,reverse_paired_merged.fq,forward_unpaired_merged.fq CC --useFwdUnpaired --FLASH --quality

contig_list.txt: unpairedDict.txt
	python dictReader.py unpairedDict.txt --runSpades --HPCC --quality
	rm -r spades_output
	rm fasta.fa

contig_lengths: contig_list.txt
	grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths
