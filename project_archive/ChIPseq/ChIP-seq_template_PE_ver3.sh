#!/bin/sh
# shellcheck disable=all

# based on mix of Yulian Tang and Soohan scripts
# the FILE var will have the full name of x.fastq.gz to remove it , uses ${FILE/.fastq.gz/} or to replace with .sam ${FILE/.fastq.gz/.sam}
# ChIP seq of CTCF active motive from BC mice CTX nuclei and HC 9/13/17
# added multiqc to last step...make changes to actual directories used.
# 12/5/19 add newer Fastqc version fastqc-0.11.8 that can handle Nova-seq better
# 11/10/21 add $HOMERVER, define $BASE, $FASTQCVER, $BOWTIEVER

BASE=~/work/ChIP/LY_MCF7_ChIP_ERdeIDR_SD_CBP_GATA3_211104
MULTIQCIN=${BASE}
MULTIQCOUT=${BASE}/analysis

DATADIR=${BASE}/data
WORKDIR=${BASE}/analysis/homer/
FQOUTDIR=${BASE}/analysis/fastqc/
BOWTIEVER=/mnt/share/apps/bio/bowtie2-2.4.5/bowtie2
HOMERVER=/opt/apps/bio/homer-4.11.1/bin
Nice 
GENOME=hg19
CPU=32

for FILE in $(ls $DATADIR/*.fastq.gz); do
echo working directory is $WORKDIR

fbase=$(basename "$FILE") 
fbase=${fbase/.fastq.gz/} 
# echo file name base is $fbase

cd ${WORKDIR}
     	nice  $FASTQCVER ${FILE} -o $FQOUTDIR
 echo    ______step0: fastqc analysis generated! ______

# run as SE first if needed  
#        nice  /opt/apps/bio/bowtie2-2.2.6/bowtie2 -q --very-sensitive -p 32 -x /mnt/share/archive/ref/align/bowtie2/$GENOME/$GENOME -S ${WORKDIR}${fbase}.sam -U ${FILE} 2> bowtie.${fbase}.log

done

# run as PE: 

# exmple of PE file names from iSeq machine
#YTC563_S1_L001_R1_001.fastq.gz YTC563_S1_L001_R2_001.fastq.gz
# or renamed files from NOva-seq:
# AG541_ChIP_CDK9cellsig_293T_191127_R2.fastq.gz

for FILE in $(ls $DATADIR/*_R1.fastq.gz); do
fbase=$(basename "$FILE") 
fbase=${fbase/_R1.fastq.gz/} 
echo "for $FILE"
echo "R1 is ${fbase}_R1.fastq.gz"
echo "R2 is ${fbase}_R2.fastq.gz"


#paired end
          nice  $BOWTIEVER -q --very-sensitive -p $CPU -x /mnt/share/archive/ref/align/bowtie2/$GENOME/$GENOME -1 ${DATADIR}/${fbase}_R1.fastq.gz -2 ${DATADIR}/${fbase}_R2.fastq.gz -S ${fbase}.sam 2> bowtie.paired.${fbase}.log


#not needed for Homer 		nice /opt/apps/bio/samtools-1.3.1/bin/samtools view -bS ${FILE/.fastq.gz/.sam} > ${FILE/.fastq.gz/_unSorted.bam}		
#not needed for Homer 		nice /opt/apps/bio/samtools-1.3.1/bin/samtools sort ${FILE/.fastq.gz/_unSorted.bam} -o ${FILE/.fastq.gz/.sorted.bam} 		
#not needed for Homer 		nice /opt/apps/bio/samtools-1.3.1/bin/samtools index ${FILE/.fastq.gz/.sorted.bam} 
#not needed for Homer		nice /opt/apps/bio/samtools-1.3.1/bin/samtools view -h ${FILE/.fastq.gz/.sorted.bam} > ${FILE/.fastq.gz/.unsorted.sam}
#not needed for Homer		nice sort -s -k 1,1 ${FILE/.fastq.gz/.unsorted.sam} > ${FILE/.fastq.gz/.sorted.sam}
# echo    ______step1: bowtie2 alignment successful ______

#not needed for Homer		nice /opt/apps/bio/homer-4.8.2/bin/makeTagDirectory ${FILE/.fastq.gz/} -format sam -unique -tbp 1 ${FILE/.fastq.gz/.sorted.sam} 
		nice $HOMERVER/makeTagDirectory $fbase -format sam -unique -tbp 1 ${WORKDIR}${fbase}.sam 
			#/opt/apps/bio/homer-4.8.2/bin/removeOutOfBoundsReads.pl ${FILE/.fastq.gz/} mm8

echo    ______step2:  makeTagDirectory of FILE.sam successful ______

     nice $HOMERVER/makeUCSCfile ${WORKDIR}${fbase} -o auto

echo    ______step3: UCSC bedgraph.gz completed ______

      nice $HOMERVER/findPeaks ${WORKDIR}${fbase} -style factor -o auto 
echo    ______step4: peak finding completed ______       
       
		nice $HOMERVER/pos2bed.pl ${WORKDIR}${fbase}/peaks.txt > ${WORKDIR}${fbase}/${fbase}.peakfile.bed
echo    ______step5: created peakfile.bed completed ______   
      
       nice $HOMERVER/annotatePeaks.pl ${WORKDIR}${fbase}/peaks.txt ${GENOME} > ${WORKDIR}${fbase}/${fbase}.peakannotation.txt 

echo    ______step6: annotation of peaks are completed ______      
       
#can use mm8 instaed of mm8r (masked) if needed.
# motif finding takes the longest, even with 12 cores. better to run searatly after finshing all other tasks..
#     nice /opt/apps/bio/homer-4.8.2/bin/findMotifsGenome.pl ${FILE/.fastq.gz/}/peaks.txt mm8 ${FILE/.fastq.gz/.MotifOutput/} -p 12 -size 200 -len 8,9,10,11,12,13,14,15 -S 200
#echo    ______step7: motifs finding completed ______ 

		nice $HOMERVER/makeBigWig.pl ${fbase} ${GENOME} -normal -force -webdir /mnt/share/archive/res/users/agamliel/ucsc/ -url  http://rosenfeldlab.ucsd.edu/res/agamliel/ucsc/

echo    ______step8: makeBigWig completed ______ 
# rm -rf ${FILE}
rm -rf ${WORKDIR}*.sam
# rm -rf ${FILE/.fastq.gz/_unSorted.bam}
# rm -rf ${FILE/.fastq.gz/.unsorted.sam}
# rm ${FILE/.fastq.gz/.sorted.bam}
# rm -rf ${FILE/.fastq.gz/.sorted.bam.bai}
# rm ${FILE/.fastq.gz/.sorted.sam}

done

 multiqc -f -v ${MULTIQCIN} -o ${MULTIQCOUT}

# Add script to collect stats from Homer libraries
#mod 2/1/18

echo "Getting Homer stats"

cd ${WORKDIR}

# Prepare headers from tmp files for Homer.run.stats.txt:

# Print the working directory on first line, next line the date:
printf   '%s\n' "Get stats from  ${WORKDIR}"  "$(date)"  >tmp1.txt

# Print tab delimited "name" "Unique Positions" "Total Tags" "Homer Peaks" "IP efficiency", and carriage return into file :
printf '%s\t%s\t%s\t%s\t%s\n' "Name" "Unique Positions" "Total Tags" "Homer Peaks" "IP efficiency" > tmp2.txt

# Make header file:
cat tmp1.txt tmp2.txt > Homer.run.stats.txt
rm tmp?.txt 


#get the number or reads used by Homer, number of peaks, IP efficiency 
for FILE in $( find ${WORKDIR}* -maxdepth 0 -type d ); do

        fbase=$(basename "$FILE") 
        # -v is variable, use to define filename.
        #NR==2 means extract line 2 only        
        awk -v n=$fbase 'NR==2 {print n,"\t",$2,"\t",$3}' ${WORKDIR}${fbase}/tagInfo.txt > tmp.txt  
        # Extract number of peaks (line 5) from peaks.txt:
        awk 'NR==5 {print $5}' ${WORKDIR}${fbase}/peaks.txt > tmp1.txt 
        # Extract  IP efficiency (line 13) from peaks.txt:
        awk 'NR==13 {print $6}' ${WORKDIR}${fbase}/peaks.txt > tmp2.txt  
        #put together so that that file tmp3.txt has 5 columns:  
        paste tmp.txt tmp1.txt tmp2.txt > tmp3.txt
        cat tmp3.txt >> Homer.run.stats.txt
        rm tmp.txt tmp?.txt
done

echo    _________ All Done ________________


