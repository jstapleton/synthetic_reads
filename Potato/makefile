# This makefile will reproduce the Potato sythetic read assembly.
# Typically, the pairedDict.txt and contig_list.txt steps are run separately as queued jobs on a computer cluster.
# qsub scripts for these two steps are provided as HPCC_Hasher.sub and HPCC_Spades.sub
# Prior to running HPCC_Spades.sub, run:
#	mkdir splitDicts
#	cd splitDicts
#	python ~/source/dictSplitter.py ../unpairedDict.txt 100
#	cd ..



all: forward_paired.fq reverse_paired.fq forward_unpaired.fq unpairedDict.txt contig_list.txt contig_lengths


clean: 
	-rm -f forward_paired.fq reverse_paired.fq forward_unpaired.fq barcodePairs.txt confirmedPairs.txt unpairedDict.txt contig_list.txt contig_lengths

forward_paired.fq reverse_paired.fq forward_unpaired.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 P_TTAGGC_L001_R1_001.fastq.gz P_TTAGGC_L001_R2_001.fastq.gz forward_paired.fq forward_unpaired.fq reverse_paired.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20
	rm reverse_unpaired.fq

unpairedDict.txt: forward_paired.fq reverse_paired.fq forward_unpaired.fq
	python barcodeHasher.py forward_paired.fq,reverse_paired.fq,forward_unpaired.fq CC --useFwdUnpaired --FLASH --quality

contig_list.txt: unpairedDict.txt
	python dictReader.py unpairedDict.txt --runSpades --HPCC --quality
	rm -r spades_output
	rm fasta.fa

contig_lengths: contig_list.txt
	grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths
