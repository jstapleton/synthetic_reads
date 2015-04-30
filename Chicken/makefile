# This makefile will reproduce the Chicken sythetic read assembly.
# Typically, the unpairedDict.txt and contig_list.txt steps are run separately as queued jobs on a computer cluster.
# qsub scripts for these two steps are provided as HPCC_Hasher.sub and HPCC_Spades.sub
# Prior to running HPCC_Spades.sub, run:
#
#	python dictSplitter.py unpairedDict.txt 4
#	for i in {1..4} ; do mkdir split$i ; mv splitDict$i.txt split$i/split$i.txt ; done
#	for i in {1..4} ; do cd split$i ; mkdir splitDicts; cd splitDicts ; python /mnt/home/jas/source/dictSplitter.py ../split$i.txt 900 ; cd ../.. ; done
#
# Then run HPCC_Spades for each of the four directories split1 - split4



all: forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq unpairedDict.txt contig_list.txt contig_lengths


clean: 
	-rm -f forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq unpairedDict.txt contig_list.txt contig_lengths barcodePairs.txt confirmedPairs.txt

forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Chicken_S1_CACTCA_L001_R1_001.fastq.gz Chicken_S1_CACTCA_L001_R2_001.fastq.gz forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Chicken_S1_CACTCA_L002_R1_001.fastq.gz Chicken_S1_CACTCA_L002_R2_001.fastq.gz forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Ch_S_E_CAGGCG_L002_R1_001.fastq.gz Ch_S_E_CAGGCG_L002_R2_001.fastq.gz forward_paired3.fq forward_unpaired3.fq reverse_paired3.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq:
	java -jar trimmomatic-0.32.jar PE -threads 20 Ch_B_E_CACTCA_L002_R1_001.fastq.gz Ch_B_E_CACTCA_L002_R2_001.fastq.gz forward_paired4.fq forward_unpaired4.fq reverse_paired4.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20
	rm reverse_unpaired.fq

forward_paired_merged.fq: forward_paired1.fq forward_paired2.fq forward_paired3.fq forward_paired4.fq
	cat forward_paired1.fq forward_paired2.fq forward_paired3.fq forward_paired4.fq > forward_paired_merged.fq

forward_unpaired_merged.fq: forward_unpaired1.fq forward_unpaired2.fq forward_unpaired3.fq forward_unpaired4.fq
	cat forward_unpaired1.fq forward_unpaired2.fq forward_unpaired3.fq forward_unpaired4.fq > forward_unpaired_merged.fq

reverse_paired_merged.fq: cat reverse_paired1.fq reverse_paired2.fq reverse_paired3.fq reverse_paired4.fq
	cat reverse_paired1.fq reverse_paired2.fq reverse_paired3.fq reverse_paired4.fq > reverse_paired_merged.fq

unpairedDict.txt: forward_paired.fq reverse_paired.fq forward_unpaired.fq
	python barcodeHasher.py forward_paired_merged.fq,reverse_paired_merged.fq,forward_unpaired_merged.fq CC --useFwdUnpaired --FLASH --quality

contig_list.txt: unpairedDict.txt
	python dictReader.py unpairedDict.txt --runSpades --HPCC --quality
	rm -r spades_output
	rm fasta.fa

contig_lengths: contig_list.txt
	grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths
