help:
	@cat Makefile

######################### WES #########################

gene_exons_pad20:
	NGSDExportGenes | cut -f1 | grep -v "#symbol" | GenesToBed -mode exon -source ensembl | BedExtend -n 20 | BedMerge -out gene_exons_pad20.bed

chunks:
	BedChunk -in gene_exons_pad20.bed -n 150 -out chunks.bed

fastq:
	BedToFasta -in chunks.bed | php fasta2fastq.php | gzip -1 > data.fastq.gz

mapped:
	/mnt/storage1/share/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 5 /tmp/local_ngs_data_GRCh38/GRCh38.fa -M data.fastq.gz | /mnt/storage1/share/opt/samtools-1.15.1/samtools view -1 - > tmp.bam
	/mnt/storage1/share/opt/samtools-1.15.1/samtools sort tmp.bam -T sort_tmp -O bam -o mapped.bam
	samtools index mapped.bam
	rm tmp.bam

mapq:
	/mnt/storage1/share/opt/samtools-1.15.1/samtools view mapped.bam | cut -f1,3,4,5 | php mapq2bed.php > mapq.bed
	TsvFilter -numeric -filter "4 = 0" -in mapq.bed  | grep -v "#" > wes_mapq_eq0.bed
	BedMerge -in wes_mapq_eq0.bed -out wes_mapq_eq0.bed


######################### WES #########################

wgs_target:
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_grch38.bed -in2 ../nobase_regions.bed -out wgs_target.bed

wgs_chunks:
	BedChunk -in wgs_target.bed -n 175 -out wgs_chunks.bed

wgs_data:
	BedToFasta -in wgs_chunks.bed | php fasta2fastq.php | gzip -1 > wgs_data.fastq.gz

wgs_mapped:
	/mnt/storage1/share/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 5 /tmp/local_ngs_data/GRCh38.fa -M wgs_data.fastq.gz | /mnt/storage1/share/opt/samtools-1.15.1/samtools view -1 - > tmp.bam
	/mnt/storage1/share/opt/samtools-1.15.1/samtools sort tmp.bam -T sort_tmp -O bam -o wgs_mapped.bam
	samtools index wgs_mapped.bam
	rm tmp.bam

wgs_mapq:
	/mnt/storage1/share/opt/samtools-1.15.1/samtools view wgs_mapped.bam | cut -f1,3,4,5 | php mapq2bed.php > wgs_mapq.bed
	TsvFilter -numeric -filter "4 = 0" -in wgs_mapq.bed  | grep -v "#" > wgs_mapq_eq0.bed
	BedMerge -in wgs_mapq_eq0.bed -out wgs_mapq_eq0.bed
