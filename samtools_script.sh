#! /bin/bash
#samtools script - using samtools for generating a VCF file for analysis

#sam to bam conversion
cd ~/temp/bowtie_output

for sample in *.sam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/.sam//')
	echo $describer
	
	#convert file from SAM to BAM format
	samtools view -bS $sample -o ${describer}.uns.bam

	#Sort BAM file
	samtools sort ${describer}.uns.bam ${describer}

	#Index BAM file
	samtools index ${describer}.bam

	#Remove intermediate files
	rm ${describer}.uns.bam

done

#getting a list of sam files
ls *.bam > Hvir_BamFiles.txt

#running samtools version 0.1.18 - getting bcf file with snp calling (-c) with potential variant sites (-v), genotypes (-g)l
samtools mpileup -uD -b Hvir_BamFiles.txt -d 70 -f ./Postjelly_capturedonly_asm_k63_05092016_FINAL_wrap.fasta -g | bcftools view -bcvg - > ./samtoolsANDvcftools_output/FieldHvir.bcf

bcftools view ./samtoolsANDvcftools_output/FieldHvir.bcf > ./samtoolsANDvcftools_output/FieldHvir.vcf #conversion to vcf format
