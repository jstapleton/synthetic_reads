


##Multiplex data analysis with two-file barcode pairing:

Trimmomatic trimming to remove low-quality sequences and adapter sequences:

```
java -jar trimmomatic-0.32.jar PE -threads 20 Lenski_F_ATCACG_L001_R1_001.fastq.gz Lenski_F_ATCACG_L001_R2_001.fastq.gz forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired1.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
```
```
java -jar trimmomatic-0.32.jar PE -threads 20 Lenski_F_ATCACG_L001_R1_001.fastq.gz Lenski_F_ATCACG_L001_R2_001.fastq.gz forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired1.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20
```
```
java -jar trimmomatic-0.32.jar PE -threads 20 Lenski_F_ATCACG_L002_R1_001.fastq.gz Lenski_F_ATCACG_L002_R2_001.fastq.gz forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired2.fq ILLUMINACLIP:JAStrim.fa:2:30:6:1:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20
```

*Output*
```
Input Read Pairs: 121901616 Both Surviving: 113753065 (93.32%) Forward Only Surviving: 5801395 (4.76%) Reverse Only Surviving: 1120970 (0.92%) Dropped: 1226186 (1.01%)
```
```
Input Read Pairs: 116841309 Both Surviving: 108843881 (93.16%) Forward Only Surviving: 5636763 (4.82%) Reverse Only Surviving: 1214564 (1.04%) Dropped: 1146101 (0.98%)
```


Demultiplex the reads:
```
mkdir demultiplexed1 demultiplexed2
cd demultiplexed1
python demultiplexer.py ../forward_paired1.fq,../reverse_paired1.fq,../forward_unpaired1.fq
cd ../demultiplexed2
python demultiplexer.py ../forward_paired2.fq,../reverse_paired2.fq,../forward_unpaired2.fq
for i in {1..24}; do cat demultiplexed1/forward_paired_"$i".fq demultiplexed2/forward_paired_"$i".fq > forward_paired"$i".fq; done
for i in {1..24}; do cat demultiplexed1/forward_unpaired_"$i".fq demultiplexed2/forward_unpaired_"$i".fq > forward_unpaired"$i".fq; done
for i in {1..24}; do cat demultiplexed1/reverse_paired_"$i".fq demultiplexed2/reverse_paired_"$i".fq > reverse_paired"$i".fq; done
```





Sort the reads according to barcode with the following HPCC .sub script:
```
qsub HPCC_Hasher.sub
```

```
#!/bin/bash -login
#PBS -l walltime=03:00:00,nodes=1:ppn=2
#PBS -j oe
#PBS -N HPC_Hasher_Flash

# Memory needed by the job
#PBS -l mem=20gb

#PBS -t 1-24
### load necessary modules, e.g.
module load Biopython

# change to the working directory where your code is located
cd $PBS_O_WORKDIR
cd strain_${PBS_ARRAYID}

# call your executable
/usr/bin/time -v python $HOME/source/barcodeHasher.py forward_paired${PBS_ARRAYID}.fq,reverse_paired${PBS_ARRAYID}.fq,forward_unpaired${PBS_ARRAYID}.fq,../Lenski_PB.fastq ATCACGC,CGATGTC,TTAGGCC,TGACCAC,ACAGTGC,GCCAATC,CAGATCC,ACTTGAC,GATCAGC,TAGCTTC,GGCTACC,CTTGTAC,AGTCAAC,AGTTCCC,ATGTCAC,CCGTCCC,GTAGAGC,GTCCGCC,GTGAAAC,GTGGCCC,GTTTCGC,CGTACGC,GAGTGGC,GGTAGCC --pairSeparateFile --useFwdUnpaired --FLASH --quality
```


```
grep -A4 'Read combination statistics' HPC_Hasher_Flash.o*
```
```
HPC_Hasher_Flash.o22607474-1:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-1-[FLASH]     Total pairs:      8498998
HPC_Hasher_Flash.o22607474-1-[FLASH]     Combined pairs:   6872512
HPC_Hasher_Flash.o22607474-1-[FLASH]     Uncombined pairs: 1626486
HPC_Hasher_Flash.o22607474-1-[FLASH]     Percent combined: 80.86%
--
HPC_Hasher_Flash.o22607474-10:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-10-[FLASH]     Total pairs:      6995189
HPC_Hasher_Flash.o22607474-10-[FLASH]     Combined pairs:   6001666
HPC_Hasher_Flash.o22607474-10-[FLASH]     Uncombined pairs: 993523
HPC_Hasher_Flash.o22607474-10-[FLASH]     Percent combined: 85.80%
--
HPC_Hasher_Flash.o22607474-11:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-11-[FLASH]     Total pairs:      7778014
HPC_Hasher_Flash.o22607474-11-[FLASH]     Combined pairs:   6465025
HPC_Hasher_Flash.o22607474-11-[FLASH]     Uncombined pairs: 1312989
HPC_Hasher_Flash.o22607474-11-[FLASH]     Percent combined: 83.12%
--
HPC_Hasher_Flash.o22607474-12:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-12-[FLASH]     Total pairs:      8514079
HPC_Hasher_Flash.o22607474-12-[FLASH]     Combined pairs:   7064693
HPC_Hasher_Flash.o22607474-12-[FLASH]     Uncombined pairs: 1449386
HPC_Hasher_Flash.o22607474-12-[FLASH]     Percent combined: 82.98%
--
HPC_Hasher_Flash.o22607474-13:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-13-[FLASH]     Total pairs:      9272393
HPC_Hasher_Flash.o22607474-13-[FLASH]     Combined pairs:   7430552
HPC_Hasher_Flash.o22607474-13-[FLASH]     Uncombined pairs: 1841841
HPC_Hasher_Flash.o22607474-13-[FLASH]     Percent combined: 80.14%
--
HPC_Hasher_Flash.o22607474-14:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-14-[FLASH]     Total pairs:      10010186
HPC_Hasher_Flash.o22607474-14-[FLASH]     Combined pairs:   8636279
HPC_Hasher_Flash.o22607474-14-[FLASH]     Uncombined pairs: 1373907
HPC_Hasher_Flash.o22607474-14-[FLASH]     Percent combined: 86.27%
--
HPC_Hasher_Flash.o22607474-15:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-15-[FLASH]     Total pairs:      13507015
HPC_Hasher_Flash.o22607474-15-[FLASH]     Combined pairs:   11022521
HPC_Hasher_Flash.o22607474-15-[FLASH]     Uncombined pairs: 2484494
HPC_Hasher_Flash.o22607474-15-[FLASH]     Percent combined: 81.61%
--
HPC_Hasher_Flash.o22607474-16:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-16-[FLASH]     Total pairs:      7702329
HPC_Hasher_Flash.o22607474-16-[FLASH]     Combined pairs:   6316390
HPC_Hasher_Flash.o22607474-16-[FLASH]     Uncombined pairs: 1385939
HPC_Hasher_Flash.o22607474-16-[FLASH]     Percent combined: 82.01%
--
HPC_Hasher_Flash.o22607474-17:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-17-[FLASH]     Total pairs:      3917164
HPC_Hasher_Flash.o22607474-17-[FLASH]     Combined pairs:   3348344
HPC_Hasher_Flash.o22607474-17-[FLASH]     Uncombined pairs: 568820
HPC_Hasher_Flash.o22607474-17-[FLASH]     Percent combined: 85.48%
--
HPC_Hasher_Flash.o22607474-18:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-18-[FLASH]     Total pairs:      8690917
HPC_Hasher_Flash.o22607474-18-[FLASH]     Combined pairs:   7248089
HPC_Hasher_Flash.o22607474-18-[FLASH]     Uncombined pairs: 1442828
HPC_Hasher_Flash.o22607474-18-[FLASH]     Percent combined: 83.40%
--
HPC_Hasher_Flash.o22607474-19:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-19-[FLASH]     Total pairs:      1982210
HPC_Hasher_Flash.o22607474-19-[FLASH]     Combined pairs:   1566517
HPC_Hasher_Flash.o22607474-19-[FLASH]     Uncombined pairs: 415693
HPC_Hasher_Flash.o22607474-19-[FLASH]     Percent combined: 79.03%
--
HPC_Hasher_Flash.o22607474-2:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-2-[FLASH]     Total pairs:      11334696
HPC_Hasher_Flash.o22607474-2-[FLASH]     Combined pairs:   9301073
HPC_Hasher_Flash.o22607474-2-[FLASH]     Uncombined pairs: 2033623
HPC_Hasher_Flash.o22607474-2-[FLASH]     Percent combined: 82.06%
--
HPC_Hasher_Flash.o22607474-20:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-20-[FLASH]     Total pairs:      6589080
HPC_Hasher_Flash.o22607474-20-[FLASH]     Combined pairs:   5476772
HPC_Hasher_Flash.o22607474-20-[FLASH]     Uncombined pairs: 1112308
HPC_Hasher_Flash.o22607474-20-[FLASH]     Percent combined: 83.12%
--
HPC_Hasher_Flash.o22607474-21:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-21-[FLASH]     Total pairs:      7466697
HPC_Hasher_Flash.o22607474-21-[FLASH]     Combined pairs:   6265378
HPC_Hasher_Flash.o22607474-21-[FLASH]     Uncombined pairs: 1201319
HPC_Hasher_Flash.o22607474-21-[FLASH]     Percent combined: 83.91%
--
HPC_Hasher_Flash.o22607474-22:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-22-[FLASH]     Total pairs:      1023717
HPC_Hasher_Flash.o22607474-22-[FLASH]     Combined pairs:   877750
HPC_Hasher_Flash.o22607474-22-[FLASH]     Uncombined pairs: 145967
HPC_Hasher_Flash.o22607474-22-[FLASH]     Percent combined: 85.74%
--
HPC_Hasher_Flash.o22607474-23:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-23-[FLASH]     Total pairs:      8283829
HPC_Hasher_Flash.o22607474-23-[FLASH]     Combined pairs:   5118968
HPC_Hasher_Flash.o22607474-23-[FLASH]     Uncombined pairs: 3164861
HPC_Hasher_Flash.o22607474-23-[FLASH]     Percent combined: 61.79%
--
HPC_Hasher_Flash.o22607474-24:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-24-[FLASH]     Total pairs:      6328953
HPC_Hasher_Flash.o22607474-24-[FLASH]     Combined pairs:   5363525
HPC_Hasher_Flash.o22607474-24-[FLASH]     Uncombined pairs: 965428
HPC_Hasher_Flash.o22607474-24-[FLASH]     Percent combined: 84.75%
--
HPC_Hasher_Flash.o22607474-3:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-3-[FLASH]     Total pairs:      10407380
HPC_Hasher_Flash.o22607474-3-[FLASH]     Combined pairs:   8727424
HPC_Hasher_Flash.o22607474-3-[FLASH]     Uncombined pairs: 1679956
HPC_Hasher_Flash.o22607474-3-[FLASH]     Percent combined: 83.86%
--
HPC_Hasher_Flash.o22607474-4:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-4-[FLASH]     Total pairs:      7524924
HPC_Hasher_Flash.o22607474-4-[FLASH]     Combined pairs:   6044779
HPC_Hasher_Flash.o22607474-4-[FLASH]     Uncombined pairs: 1480145
HPC_Hasher_Flash.o22607474-4-[FLASH]     Percent combined: 80.33%
--
HPC_Hasher_Flash.o22607474-5:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-5-[FLASH]     Total pairs:      10519817
HPC_Hasher_Flash.o22607474-5-[FLASH]     Combined pairs:   8915165
HPC_Hasher_Flash.o22607474-5-[FLASH]     Uncombined pairs: 1604652
HPC_Hasher_Flash.o22607474-5-[FLASH]     Percent combined: 84.75%
--
HPC_Hasher_Flash.o22607474-6:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-6-[FLASH]     Total pairs:      6322361
HPC_Hasher_Flash.o22607474-6-[FLASH]     Combined pairs:   5117747
HPC_Hasher_Flash.o22607474-6-[FLASH]     Uncombined pairs: 1204614
HPC_Hasher_Flash.o22607474-6-[FLASH]     Percent combined: 80.95%
--
HPC_Hasher_Flash.o22607474-7:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-7-[FLASH]     Total pairs:      13138610
HPC_Hasher_Flash.o22607474-7-[FLASH]     Combined pairs:   10980757
HPC_Hasher_Flash.o22607474-7-[FLASH]     Uncombined pairs: 2157853
HPC_Hasher_Flash.o22607474-7-[FLASH]     Percent combined: 83.58%
--
HPC_Hasher_Flash.o22607474-8:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-8-[FLASH]     Total pairs:      8453463
HPC_Hasher_Flash.o22607474-8-[FLASH]     Combined pairs:   7076962
HPC_Hasher_Flash.o22607474-8-[FLASH]     Uncombined pairs: 1376501
HPC_Hasher_Flash.o22607474-8-[FLASH]     Percent combined: 83.72%
--
HPC_Hasher_Flash.o22607474-9:[FLASH] Read combination statistics:
HPC_Hasher_Flash.o22607474-9-[FLASH]     Total pairs:      8935049
HPC_Hasher_Flash.o22607474-9-[FLASH]     Combined pairs:   7645110
HPC_Hasher_Flash.o22607474-9-[FLASH]     Uncombined pairs: 1289939
HPC_Hasher_Flash.o22607474-9-[FLASH]     Percent combined: 85.56%
```



```
grep -A1 'total read pairs' HPC_Hasher_Flash.o*
```
*Returned stats:*
```
HPC_Hasher_Flash.o22607474-1:8939290 total read pairs, 8872508 passed, 66782 non-compliant barcodes, AB counts: [8872508, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-1-Found 529 barcode pairs
--
HPC_Hasher_Flash.o22607474-10:7274242 total read pairs, 7254760 passed, 19482 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 1, 0, 0, 7254759, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-10-Found 57 barcode pairs
--
HPC_Hasher_Flash.o22607474-11:8151595 total read pairs, 8048677 passed, 102918 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8048677, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-11-Found 921 barcode pairs
--
HPC_Hasher_Flash.o22607474-12:8979352 total read pairs, 8929129 passed, 50223 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8929128, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-12-Found 1509 barcode pairs
--
HPC_Hasher_Flash.o22607474-13:9838249 total read pairs, 9800591 passed, 37658 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 9800586, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-13-Found 915 barcode pairs
--
HPC_Hasher_Flash.o22607474-14:10373402 total read pairs, 10318263 passed, 55139 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 10318245, 0, 1, 0, 0, 0, 0, 16, 0, 0, 0]
HPC_Hasher_Flash.o22607474-14-Found 672 barcode pairs
--
HPC_Hasher_Flash.o22607474-15:14377688 total read pairs, 14214396 passed, 163292 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 14214394, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-15-Found 1440 barcode pairs
--
HPC_Hasher_Flash.o22607474-16:8106004 total read pairs, 8038491 passed, 67513 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8038489, 0, 0, 0, 0, 0, 1, 0, 0]
HPC_Hasher_Flash.o22607474-16-Found 495 barcode pairs
--
HPC_Hasher_Flash.o22607474-17:4089676 total read pairs, 4071121 passed, 18555 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4071121, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-17-Found 197 barcode pairs
--
HPC_Hasher_Flash.o22607474-18:9095791 total read pairs, 9047392 passed, 48399 non-compliant barcodes, AB counts: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9047391, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-18-Found 1023 barcode pairs
--
HPC_Hasher_Flash.o22607474-19:2105827 total read pairs, 2093789 passed, 12038 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2093789, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-19-Found 261 barcode pairs
--
HPC_Hasher_Flash.o22607474-2:11992485 total read pairs, 11930637 passed, 61848 non-compliant barcodes, AB counts: [0, 11930637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-2-Found 662 barcode pairs
--
HPC_Hasher_Flash.o22607474-20:6892344 total read pairs, 6834299 passed, 58045 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 6834295, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-20-Found 647 barcode pairs
--
HPC_Hasher_Flash.o22607474-21:7847884 total read pairs, 7798683 passed, 49201 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 7798681, 0, 0, 0]
HPC_Hasher_Flash.o22607474-21-Found 661 barcode pairs
--
HPC_Hasher_Flash.o22607474-22:1067098 total read pairs, 1060475 passed, 6623 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1060475, 0, 0]
HPC_Hasher_Flash.o22607474-22-Found 136 barcode pairs
--
HPC_Hasher_Flash.o22607474-23:8724173 total read pairs, 8673609 passed, 50564 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8673608, 0]
HPC_Hasher_Flash.o22607474-23-Found 477 barcode pairs
--
HPC_Hasher_Flash.o22607474-24:6625000 total read pairs, 6553945 passed, 71055 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6553945]
HPC_Hasher_Flash.o22607474-24-Found 834 barcode pairs
--
HPC_Hasher_Flash.o22607474-3:10876544 total read pairs, 10809493 passed, 67051 non-compliant barcodes, AB counts: [0, 0, 10809493, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-3-Found 1532 barcode pairs
--
HPC_Hasher_Flash.o22607474-4:7996412 total read pairs, 7954228 passed, 42184 non-compliant barcodes, AB counts: [0, 0, 0, 7954228, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-4-Found 831 barcode pairs
--
HPC_Hasher_Flash.o22607474-5:10976255 total read pairs, 10911001 passed, 65254 non-compliant barcodes, AB counts: [0, 0, 0, 0, 10911001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-5-Found 1674 barcode pairs
--
HPC_Hasher_Flash.o22607474-6:6701926 total read pairs, 6665185 passed, 36741 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 6665184, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-6-Found 380 barcode pairs
--
HPC_Hasher_Flash.o22607474-7:13778540 total read pairs, 13681197 passed, 97343 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 13681197, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-7-Found 1539 barcode pairs
--
HPC_Hasher_Flash.o22607474-8:8889354 total read pairs, 8853848 passed, 35506 non-compliant barcodes, AB counts: [0, 0, 0, 0, 0, 0, 0, 8853848, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-8-Found 1104 barcode pairs
--
HPC_Hasher_Flash.o22607474-9:9356523 total read pairs, 9302048 passed, 54475 non-compliant barcodes, AB counts: [0, 0, 0, 1, 0, 0, 0, 0, 9302047, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
HPC_Hasher_Flash.o22607474-9-Found 125 barcode pairs
```



Break the resulting JSON-dumped dictionary into smaller dictionaries to allow parallel assembly on a cluster:

```
mkdir strain_{1..24}/splitDicts
for i in {1..24}; do cd strain_$i/splitDicts ; python ~/source/dictSplitter.py ../pairedDict.txt 200 ; cd ../.. ; done
```


Assemble in parallel on a cluster with the following submission script 
(works on MSU's HPCC, may need to be modified for other clusters):

Submit with:
```
for i in {1..24} ; do cd strain_$i ; qsub ../HPCC_Spades.sub ; cd .. ; done
for i in {1..4} ; do cd strain_$i ; qsub ../HPCC_Spades.sub ; cd .. ; done
for i in {5..8} ; do cd strain_$i ; qsub ../HPCC_Spades.sub ; cd .. ; done
for i in {14..18} ; do cd strain_$i ; qsub ../HPCC_Spades.sub ; cd .. ; done

```





```
#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=2
#PBS -j oe
#PBS -N HCT_Spades_parallel

# Memory needed by the job
#PBS -l mem=2gb

#PBS -t 1-900

### load necessary modules
module load Biopython
module swap GNU GNU/4.7.1
module load SPAdes/3.1.1

cd $PBS_O_WORKDIR
export WORKDIR="/mnt/scratch/${USER}/${PBS_JOBNAME}/"
mkdir -p $WORKDIR
cp $PBS_O_WORKDIR/splitDicts/splitDict${PBS_ARRAYID}.txt $WORKDIR
cd $WORKDIR

##do the work here
/usr/bin/time -v python $HOME/source/dictReader.py splitDict${PBS_ARRAYID}.txt --runSpades --HPCC --quality

##copy data back
cp -R $WORKDIR ${PBS_O_WORKDIR}/${PBS_JOBID}
rm -r $WORKDIR
```

Combine the resulting contig list files into one:

```
cat XXXXX[*/contig_list.txt > contig_list.txt
```
... where XXXXX is the PBS_JOBID


Pull out the lengths of the assembled synthetic long reads:

```
for i in {1..24} ; do cd strain_$i ; cat */contig_list.txt > contig_list.txt ; grep 'NODE' contig_list.txt | cut -d_ -f 4 | sort -n -r > contig_lengths ; awk '{s+=$1} END {print s}' contig_lengths ; cd ../ ; done
```


