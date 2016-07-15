# This makefile will reproduce the env sythetic read assembly.
# Typically, the unpairedDict.txt and contig_list.txt steps are run separately as queued jobs on a computer cluster.
# qsub scripts for these two steps are provided as HPCC_Hasher.sub and HPCC_Spades.sub
# Prior to running HPCC_Spades.sub, run:
#
#	mkdir splitDicts
#	cd splitDicts
#	python ~/source/dictSplitter.py ../pairedDict.txt 900
#	cd ..


all: forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq pb_reads.fq pairedDict.txt contig_list.txt contig_lengths longest_contigs.txt unique_contig_list.txt unique_contig_lengths filtered_contigs.txt filtered_contigs_noBarcodes.txt env1.water env2.water env1RC.water env2RC.water env1waterIdentity env1RCwaterIdentity env2waterIdentity env2RCwaterIdentity env1wRate env2wRate env1RCwRate env2RCwRate env1maxwRate env2maxwRate


clean: 
	-rm -f forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq forward_paired_merged.fq reverse_paired_merged.fq forward_unpaired_merged.fq pairedDict.txt contig_list.txt contig_lengths barcodePairs.txt confirmedPairs.txt longest_contigs.txt unique_contig_list.txt unique_contig_lengths filtered_contigs.txt filtered_contigs_noBarcodes.txt env1.water env2.water env1RC.water env2RC.water env1waterIdentity env1RCwaterIdentity env2waterIdentity env2RCwaterIdentity env1wRate env2wRate env1RCwRate env2RCwRate env1maxwRate env2maxwRate

forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq:
	java -jar trimmomatic-0.32.jar PE Env-F_S2_L001_R1_001.fastq.gz Env-PB_S2_L001_R2_001.fastq.gz forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq:
	java -jar trimmomatic-0.32.jar PE Env-F_S2_L002_R1_001.fastq.gz Env-F_S2_L002_R2_001.fastq.gz forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

forward_paired3.fq reverse_paired3.fq forward_unpaired3.fq:
	java -jar trimmomatic-0.32.jar PE Env-F_S2_L003_R1_001.fastq.gz Env-F_S2_L003_R2_001.fastq.gz forward_paired3.fq forward_unpaired3.fq reverse_paired3.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

forward_paired4.fq reverse_paired4.fq forward_unpaired4.fq:
	java -jar trimmomatic-0.32.jar PE Env-F_S2_L004_R1_001.fastq.gz Env-F_S2_L004_R2_001.fastq.gz forward_paired4.fq forward_unpaired4.fq reverse_paired4.fq reverse_unpaired.fq ILLUMINACLIP:JAStrim.fa:2:40:6:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
	rm reverse_unpaired.fq

forward_paired_merged.fq: forward_paired1.fq forward_paired2.fq forward_paired3.fq forward_paired4.fq
	cat forward_paired1.fq forward_paired2.fq forward_paired3.fq forward_paired4.fq > forward_paired_merged.fq

forward_unpaired_merged.fq: forward_unpaired1.fq forward_unpaired2.fq forward_unpaired3.fq forward_unpaired4.fq
	cat forward_unpaired1.fq forward_unpaired2.fq forward_unpaired3.fq forward_unpaired4.fq > forward_unpaired_merged.fq

reverse_paired_merged.fq: cat reverse_paired1.fq reverse_paired2.fq reverse_paired3.fq reverse_paired4.fq
	cat reverse_paired1.fq reverse_paired2.fq reverse_paired3.fq reverse_paired4.fq > reverse_paired_merged.fq

pb_reads.fq: Env-PB_S1_L001_R1_001.fastq Env-PB_S1_L002_R1_001.fastq Env-PB_S1_L003_R1_001.fastq Env-PB_S1_L004_R1_001.fastq
	cat Env-PB_S1_L001_R1_001.fastq Env-PB_S1_L002_R1_001.fastq Env-PB_S1_L003_R1_001.fastq Env-PB_S1_L004_R1_001.fastq > pb_reads.fq

pairedDict.txt: forward_paired.fq reverse_paired.fq forward_unpaired.fq
	python barcodeHasher.py forward_paired_merged.fq,reverse_paired_merged.fq,forward_unpaired_merged.fq,pb_reads.fq AC --useFwdUnpaired --pairSeparateFile --FLASH --quality --SCA2

contig_list.txt: pairedDict.txt
	python dictReader.py pairedDict.txt --runSpades --HPCC --quality
	rm -r spades_output
	rm fasta.fa

contig_lengths: contig_list.txt
	grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths

longest_contigs.txt: contig_list.txt
	python longest_contig_writer.py contig_list.txt

unique_contig_list: longest_contigs.txt
	python unique_barcodes.py longest_contigs.txt

unique_contig_lengths: unique_contig_list.txt
	grep 'NODE' unique_contig_list.txt | cut -d_ -f 4 | sort -n -r > unique_contig_lengths

filtered_contigs.txt: unique_contig_list.txt
	python ~/source/contigFilter.py unique_contig_list.txt --minlength 500 --maxlength 3200

filtered_contigs_noBarcodes.txt: filtered_contigs.txt
	sed '/Barcode/d' filtered_contigs.txt > filtered_contigs_noBarcodes.txt

env1.water: Env1.fa filtered_contigs_noBarcodes.txt 
	water Env1.fa filtered_contigs_noBarcodes.txt -gapopen 10.0 -gapextend 0.5 -outfile env1.water

env2.water: Env2.fa filtered_contigs_noBarcodes.txt
	water Env2.fa filtered_contigs_noBarcodes.txt -gapopen 10.0 -gapextend 0.5 -outfile env2.water

env1RC.water: Env1.fa filtered_contigs_noBarcodes.txt
	water Env1.fa filtered_contigs_noBarcodes.txt -gapopen 10.0 -gapextend 0.5 -sreverse1 -outfile env1RC.water

env2RC.water: Env2.fa filtered_contigs_noBarcodes.txt
	water Env2.fa filtered_contigs_noBarcodes.txt -gapopen 10.0 -gapextend 0.5 -sreverse1 -outfile env2RC.water

env1waterIdentity: env1.water
	grep 'Identity' env1.water | cut -d' ' -f6 > env1waterIdentity

env1RCwaterIdentity: env1RC.water
	grep 'Identity' env1RC.water | cut -d' ' -f6 > env1RCwaterIdentity

env2waterIdentity: env2.water
	grep 'Identity' env2.water | cut -d' ' -f6 > env2waterIdentity

env2RCwaterIdentity: env2RC.water
	grep 'Identity' env2RC.water | cut -d' ' -f6 > env2RCwaterIdentity

env1wRate: env1waterIdentity
	python EnvErrorRates.py env1waterIdentity > env1wRate

env2wRate: env2waterIdentity
	python EnvErrorRates.py env2waterIdentity > env2wRate

env1RCwRate: env1RCwaterIdentity
	python EnvErrorRates.py env1RCwaterIdentity > env1RCwRate

env2RCwRate: env2RCwaterIdentity
	python EnvErrorRates.py env2RCwaterIdentity > env2RCwRate

env1maxwRate: env1wRate env1RCwRate
	python listZipper.py env1wRate env1RCwRate > env1maxwRate

env2maxwRate: env2wRate env2RCwRate
	python listZipper.py env2wRate env2RCwRate > env2maxwRate
