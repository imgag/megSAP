help:
	@cat Makefile
	
vcf2gsvar:
	cp an_vep_out1.vcf vcf2gsvar_in1.vcf

vc_freebayes:
	#creating input BAM
	cut -f 1,2,3 /mnt/projects/test/KontrollDNACoriell/Sample_NA12878_13/NA12878_13.GSvar | egrep "chr3|chr7|chr8|chr10" | awk '{print $$1 "\t" $$2 "\t" ($$3+1)}' > tmp_vars.bed
	samtools view -L tmp_vars.bed /mnt/projects/test/KontrollDNACoriell/Sample_NA12878_13/NA12878_13.bam | cut -f1 | sort --uniq > tmp_reads.txt
	FastqExtract -in /mnt/projects/test/KontrollDNACoriell/Sample_NA12878_13/NA12878_13_S21_L001_R1_001.fastq.gz -ids tmp_reads.txt -out tmp_L001_R1_001.fastq.gz
	FastqExtract -in /mnt/projects/test/KontrollDNACoriell/Sample_NA12878_13/NA12878_13_S21_L001_R2_001.fastq.gz -ids tmp_reads.txt -out tmp_L001_R2_001.fastq.gz
	php /mnt/users/all/megSAP/src/NGS/mapping_bwa.php -in1 tmp_L001_R1_001.fastq.gz -in2 tmp_L001_R2_001.fastq.gz -out vc_freebayes_in.bam
	samtools index vc_freebayes_in.bam
	# creating BED file
	BedExtend -in /mnt/share/data/enrichment/ssTruRiskv1.2_2015_05_18.bed -n 50 -out tmp_roi_ext.bed
	BedIntersect -in tmp_roi_ext.bed -in2 tmp_vars.bed -mode in -out tmp_roi_overlap.bed
	BedShrink -in tmp_roi_overlap.bed -n 50 -out vc_freebayes_in.bed
	rm tmp*