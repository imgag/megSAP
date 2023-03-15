create_custom_somatic_panel:
	BedIntersect -in /mnt/storage2/megSAP/data/enrichment/ssSC_v5.bed -in2 /mnt/storage2/megSAP/data/enrichment/twistCustomExomeV2_2021_12_14.bed -mode in -out somatic_custom_panel.bed
	BedAnnotateGenes -in somatic_custom_panel.bed | cut -f4 | sort | uniq > somatic_custom_panel_genes.txt