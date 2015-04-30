# This makefile will reproduce the MG1655 sythetic read assembly.
# Typically, the pairedDict.txt and contig_list.txt steps are run separately as queued jobs on a computer cluster.
# qsub scripts for these two steps are provided as HPCC_Hasher.sub and HPCC_Spades.sub
# Prior to running HPCC_Spades.sub, run:
#	mkdir splitDicts
#	cd splitDicts
#	python ~/source/dictSplitter.py ../pairedDict.txt 200
#	cd ..



all: forward_paired.fq reverse_paired.fq forward_unpaired.fq foundPairs.fastq trimmed.fastq pairedDict.txt contig_list.txt contig_lengths MG1655.fa.amb MG1655.fa.bwt MG1655.fa.sa MG1655.fa.ann MG1655.fa.pac contigs_noBarcodes.txt alignMG1655.sam


clean: 
	-rm -f forward_paired.fq reverse_paired.fq forward_unpaired.fq foundPairs.fastq trimmed.fastq barcodePairs.txt confirmedPairs.txt pairedDict.txt contig_list.txt contig_lengths

forward_paired.fq reverse_paired.fq forward_unpaired.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 MG1655_F_ATCACG_L002_R1_001.fastq.gz MG1655_F_ATCACG_L002_R2_001.fastq.gz.fastq.gz testR.fastq.gz forward_paired.fq forward_unpaired.fq reverse_paired.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
	rm reverse_unpaired.fq

foundPairs.fastq:
	python singleFilePairing.py testF.fastq AC,CC

trimmed.fastq: foundPairs.fastq
	java -jar trimmomatic-0.32.jar SE foundPairs.fastq trimmed.fastq ILLUMINACLIP:pairedTrim.txt:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

pairedDict.txt: forward_paired.fq reverse_paired.fq forward_unpaired.fq trimmed.fastq
	python barcodeHasher.py forward_paired.fq,reverse_paired.fq,forward_unpaired.fq,trimmed.fastq CC,AC --useFwdUnpaired --pairSameFile

contig_list.txt: pairedDict.txt
	python dictReader.py pairedDict.txt --runSpades --HPCC
	rm -r spades_output
	rm fasta.fa

contig_lengths: contig_list.txt
	grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths

# Make MG1655 index for BWA-MEM alignment
MG1655.fa.amb MG1655.fa.bwt MG1655.fa.sa MG1655.fa.ann MG1655.fa.pac: MG1655.fa
	bwa index MG1655.fa

contigs_noBarcodes.txt: contig_list.txt
	sed '/Barcode/d' contig_list.txt > contigs_noBarcodes.txt

alignMG1655.sam: contigs_noBarcodes.txt
	bwa mem MG1655.fa contigs_noBarcodes.txt > alignMG1655.sam
