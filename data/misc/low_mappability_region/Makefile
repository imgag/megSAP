help:
	@cat Makefile

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
	TsvFilter -numeric -filter "4 = 0" -in mapq.bed  | grep -v "#" > mapq_eq0.bed
	BedMerge -in mapq_eq0.bed -out mapq_eq0.bed
