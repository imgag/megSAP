help:
	@cat Makefile
	
vcf2gsvar:
	cp an_vep_out1.vcf vcf2gsvar_in1.vcf
	cp an_vep_out2.vcf vcf2gsvar_in2.vcf
	cp an_vep_out_db2.vcf vcf2gsvar_in_db2.vcf
	cp an_vep_out_empty.vcf vcf2gsvar_in_empty.vcf
	cp an_vep_out_dragen.vcf vcf2gsvar_in_dragen.vcf

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

#Creating input BAM files for somatic pipeline Manta SV positive control
sv_somatic:
	samtools view -h /mnt/projects/diagnostic/SomaticAndTreatment/Sample_DX197143_01/DX197143_01.bam chr10:123238830-123240576 chr12:80078296-80080042 | sed "s/DX197143_01/manta_sv_tumor/g" | samtools view -bSh > temp_tumor.bam
	samtools view -h /mnt/projects/diagnostic/SomaticAndTreatment/Sample_DX196743_01/DX196743_01.bam chr10:123238830-123240576 chr12:80078296-80080042 | sed "s/DX196743_01/manta_sv_normal/g" | samtools view -bSh > temp_normal.bam
	samtools view -H temp_tumor.bam | sed "s/DX197143_01/manta_sv_tumor/g" > temp_header_tumor.sam
	samtools view -H temp_normal.bam | sed "s/DX196743_01/manta_sv_normal/g" > temp_header_normal.sam
	samtools reheader temp_header_tumor.sam  temp_tumor.bam > manta_sv_tumor.bam
	samtools reheader temp_header_normal.sam  temp_normal.bam > manta_sv_normal.bam 
	samtools index manta_sv_tumor.bam
	samtools index manta_sv_normal.bam
	rm temp_tumor.bam temp_header_tumor.sam temp_normal.bam temp_header_normal.sam

#Creating RNA test files for somatic pipeline test
an_somatic_gsvar_rna:
	cat /mnt/projects/diagnostic/SomaticAndTreatment/Sample_RX196019_01/RX196019_01_counts.tsv | egrep "#|TPTE|BAGE2|KRTAP13-2|CFAP298|AP000275.65|C21orf58" > somatic_pipeline_rna_counts.tsv
	samtools view -h /mnt/projects/diagnostic/SomaticAndTreatment/Sample_RX196019_01/RX196019_01.bam chr21:10971403-10971503 chr21:11039344-11039444 chr21:31744189-31744289 chr21:33974044-33974144 chr21:47737800-47737900 | sed "s/RX196019_01/somatic_rna/g" | samtools view -bSh > temp_rna.bam
	samtools view -H temp_rna.bam | sed "s/RX196019_01/somatic_rna/g" > temp_header_rna.sam
	samtools reheader temp_header_rna.sam temp_rna.bam > somatic_pipeline_tumor_rna.bam
	samtools index somatic_pipeline_tumor_rna.bam
	rm temp_rna.bam temp_header_rna.sam
vc_varscan2:
	samtools view -h -L vc_varscan2.bed /mnt/projects/diagnostic/SomaticAndTreatment/Sample_DX174130_01/DX174130_01.bam | sed "s/DX174130_01/varscan2_tumor/g" |samtools view -bSh >  vc_varscan2_tu_temp.bam
	samtools view -H vc_varscan2_tu_temp.bam | sed "s/DX174130_01/varscan2_tumor/g" > varscan_temp_header.sam
	samtools reheader varscan_temp_header.sam vc_varscan2_tu_temp.bam > vc_varscan2_tu.bam
	samtools index vc_varscan2_tu.bam
	rm varscan_temp_header.sam vc_varscan2_tu_temp.bam

an_somatic_cancerhotspots:
	zcat /mnt/projects/diagnostic/SomaticAndTreatment/Somatic_DX205275_01-DX200281_01/DX205275_01-DX200281_01_var_annotated.vcf.gz | egrep "#|NRAS|BRAF|CDKN2A|TERT|JAK2" | sed "s/DX205275_01/DX000000_01/g" | sed "s/DX200281_01/DX000000_01/g" > an_somatic_cancerhotspots_in1.vcf