help:
	@cat Makefile
	
vc_freebayes:
	#creating input BAM
	grep chr22 /mnt/storage3/projects/test/KontrollDNACoriell/Sample_NA12878_58/NA12878_58.GSvar | head -50 | cut -f 1,2,3 | awk '{print $$1 "\t" $$2 "\t" ($$3+1)}' > tmp_vars.bed
	samtools view -L tmp_vars.bed /mnt/storage3/projects/test/KontrollDNACoriell/Sample_NA12878_58/NA12878_58.bam | cut -f1 | sort --uniq > tmp_reads.txt
	FastqExtract -in /mnt/storage3/projects/test/KontrollDNACoriell/Sample_NA12878_58/NA12878_58_BamToFastq_R1_001.fastq.gz -ids tmp_reads.txt -out tmp_L001_R1_001.fastq.gz
	FastqExtract -in /mnt/storage3/projects/test/KontrollDNACoriell/Sample_NA12878_58/NA12878_58_BamToFastq_R2_001.fastq.gz -ids tmp_reads.txt -out tmp_L001_R2_001.fastq.gz
	php /mnt/storage2/megSAP/pipeline/src/NGS/mapping_bwa.php -in1 tmp_L001_R1_001.fastq.gz -in2 tmp_L001_R2_001.fastq.gz -out vc_freebayes_in.bam
	samtools index vc_freebayes_in.bam
	# creating BED file
	BedExtend -in /mnt/storage2/megSAP/data/enrichment/ssHAEv7_2021_06_08.bed -n 50 -out tmp_roi_ext.bed
	BedIntersect -in tmp_roi_ext.bed -in2 tmp_vars.bed -mode in -out tmp_roi_overlap.bed
	BedShrink -in tmp_roi_overlap.bed -n 50 -out vc_freebayes_in.bed
	rm tmp*
	
vc_mosaic:
	@cd ../ && php vc_mosaic_create_test_data.php

vcf2gsvar:
	cp an_vep_out1.vcf vcf2gsvar_in1.vcf
	cp an_vep_out_empty.vcf vcf2gsvar_in_empty.vcf
	cp an_vep_out_dragen.vcf vcf2gsvar_in_dragen.vcf
	zcat ../../test/data_multi/multi_var_annotated.vcf.gz > vcf2gsvar_in_multi.vcf

#Creating input BAM files for somatic pipeline Manta SV positive control
sv_somatic:
	samtools view -h /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX197143_01/DX197143_01.bam chr10:121479316-121481062 chr12:79684516-79686262 | sed "s/DX197143_01/manta_sv_tumor/g" | samtools view -bSh > temp_tumor.bam
	samtools view -h /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX196743_01/DX196743_01.bam chr10:121479316-121481062 chr12:79684516-79686262 | sed "s/DX196743_01/manta_sv_normal/g" | samtools view -bSh > temp_normal.bam
	samtools view -H temp_tumor.bam | sed "s/DX197143_01/manta_sv_tumor/g" > temp_header_tumor.sam
	samtools view -H temp_normal.bam | sed "s/DX196743_01/manta_sv_normal/g" > temp_header_normal.sam
	samtools reheader temp_header_tumor.sam  temp_tumor.bam > manta_sv_tumor.bam
	samtools reheader temp_header_normal.sam  temp_normal.bam > manta_sv_normal.bam 
	samtools index manta_sv_tumor.bam
	samtools index manta_sv_normal.bam
	rm temp_tumor.bam temp_header_tumor.sam temp_normal.bam temp_header_normal.sam

#Creating RNA test files for somatic pipeline test
an_somatic_gsvar_rna:
	cat /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_RX196019_01/RX196019_01_counts.tsv | egrep "#|TPTE|BAGE2|KRTAP13-2|CFAP298|AP000275.65|C21orf58" > somatic_pipeline_rna_counts.tsv
	samtools view -h /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_RX196019_01/RX196019_01.bam chr21:10473013-10473113 chr21:10540954-10541054 chr21:30371871-30371971 chr21:32601734-32601834 chr21:46317886-46317986 | sed "s/RX196019_01/somatic_rna/g" | samtools view -bSh > temp_rna.bam
	samtools view -H temp_rna.bam | sed "s/RX196019_01/somatic_rna/g" > temp_header_rna.sam
	samtools reheader temp_header_rna.sam temp_rna.bam > somatic_pipeline_tumor_rna.bam
	samtools index somatic_pipeline_tumor_rna.bam
	rm temp_rna.bam temp_header_rna.sam
vc_varscan2:
	samtools view -h -L vc_varscan2.bed /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX174130_01/DX174130_01.bam | sed "s/DX174130_01/varscan2_tumor/g" |samtools view -bSh >  vc_varscan2_tu_temp.bam
	samtools view -H vc_varscan2_tu_temp.bam | sed "s/DX174130_01/varscan2_tumor/g" > varscan_temp_header.sam
	samtools reheader varscan_temp_header.sam vc_varscan2_tu_temp.bam > vc_varscan2_tu.bam
	samtools index vc_varscan2_tu.bam
	rm varscan_temp_header.sam vc_varscan2_tu_temp.bam

an_somatic_cancerhotspots:
	zcat /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Somatic_DX205275_01-DX200281_01/DX205275_01-DX200281_01_var_annotated.vcf.gz | egrep "#|NRAS|BRAF|CDKN2A|TERT|JAK2" | sed "s/DX205275_01/DX000000_01/g" | sed "s/DX200281_01/DX000000_01/g" > an_somatic_cancerhotspots_in1.vcf
	
	
vc_strelka:
	#tumor
	samtools view -hb  /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX170388_01/DX170388_01.bam chr12:48608347-54604173 > vc_strelka2_tu_in.bam 
	samtools index vc_strelka2_tu_in.bam
	#normal
	samtools view -hb /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX170387_01/DX170387_01.bam chr12:48608347-54604173 > vc_strelka2_no_in.bam
	samtools index vc_strelka2_no_in.bam
baf_somatic:
	cp vc_strelka2_2_out.vcf.gz baf_somatic_sites.vcf.gz

hla_genotyper:
	samtools view -hb /mnt/storage3/projects/external/dante/Sample_EX000024_02/EX000024_02.bam chr6:29910331-29910403 chr6:29910534-29910803 chr6:29911045-29911320 chr6:29911899-29912174 | samtools sort > hla_genotyper_in.bam
	samtools index hla_genotyper_in.bam