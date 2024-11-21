CURDATE=$(shell date +"%Y-%m-%d")
NGSBITS=PATH_TO_NGS_BITS_TOOLS

help:
	cat Makefile

export: export_somatic export_germline export_germline_sv

export_germline:
	mkdir -p $(CURDATE)
	TMPDIR=/var/lib/mysql/tmp/ $(NGSBITS)/NGSDExportAnnotationData -germline $(CURDATE)/NGSD_germline_unsorted.vcf -genes $(CURDATE)/NGSD_genes.bed > $(CURDATE)/NGSD_germline.log 2>&1
	$(NGSBITS)/VcfStreamSort -in $(CURDATE)/NGSD_germline_unsorted.vcf -out $(CURDATE)/NGSD_germline.vcf
	bgzip -c $(CURDATE)/NGSD_germline.vcf > $(CURDATE)/NGSD_germline.vcf.gz
	tabix -p vcf $(CURDATE)/NGSD_germline.vcf.gz
	#cleanup
	rm $(CURDATE)/NGSD_germline_unsorted.vcf
	rm $(CURDATE)/NGSD_germline.vcf
	#linking
	rm -f NGSD_genes.bed
	ln -s $(CURDATE)/NGSD_genes.bed NGSD_genes.bed
	rm -f NGSD_germline.vcf.gz
	rm -f NGSD_germline.vcf.gz.tbi
	ln -s $(CURDATE)/NGSD_germline.vcf.gz NGSD_germline.vcf.gz
	ln -s $(CURDATE)/NGSD_germline.vcf.gz.tbi NGSD_germline.vcf.gz.tbi
	
export_germline_sv:
	mkdir -p $(CURDATE)
	$(NGSBITS)/NGSDExportSV -out_folder $(CURDATE) > $(CURDATE)/NGSD_germline_sv.log 2>&1
	
	#sort translocations
	$(NGSBITS)/BedpeSort -in $(CURDATE)/sv_translocation.bedpe -out $(CURDATE)/sv_translocation.bedpe
	
	#compress files
	bgzip -c $(CURDATE)/sv_deletion.bedpe > $(CURDATE)/sv_deletion.bedpe.gz
	bgzip -c $(CURDATE)/sv_duplication.bedpe > $(CURDATE)/sv_duplication.bedpe.gz
	bgzip -c $(CURDATE)/sv_insertion.bedpe > $(CURDATE)/sv_insertion.bedpe.gz
	bgzip -c $(CURDATE)/sv_inversion.bedpe > $(CURDATE)/sv_inversion.bedpe.gz
	bgzip -c $(CURDATE)/sv_translocation.bedpe > $(CURDATE)/sv_translocation.bedpe.gz
	
	#create index
	tabix -0 -b 2 -e 5 $(CURDATE)/sv_deletion.bedpe.gz
	tabix -0 -b 2 -e 5 $(CURDATE)/sv_duplication.bedpe.gz
	tabix -0 -b 2 -e 3 $(CURDATE)/sv_insertion.bedpe.gz
	tabix -0 -b 2 -e 5 $(CURDATE)/sv_inversion.bedpe.gz
	tabix -0 -b 2 -e 3 $(CURDATE)/sv_translocation.bedpe.gz
	
	#cleanup
	rm $(CURDATE)/sv_*.bedpe
	
	rm -f sv_deletion.bedpe.gz
	rm -f sv_duplication.bedpe.gz
	rm -f sv_insertion.bedpe.gz
	rm -f sv_inversion.bedpe.gz
	rm -f sv_translocation.bedpe.gz
	rm -f sv_breakpoint_density.igv
	rm -f sv_breakpoint_density_*.igv
	
	rm -f sv_deletion.bedpe.gz.tbi
	rm -f sv_duplication.bedpe.gz.tbi
	rm -f sv_insertion.bedpe.gz.tbi
	rm -f sv_inversion.bedpe.gz.tbi
	rm -f sv_translocation.bedpe.gz.tbi
	
	#link current export
	ln -s $(CURDATE)/sv_deletion.bedpe.gz sv_deletion.bedpe.gz
	ln -s $(CURDATE)/sv_duplication.bedpe.gz sv_duplication.bedpe.gz
	ln -s $(CURDATE)/sv_insertion.bedpe.gz sv_insertion.bedpe.gz
	ln -s $(CURDATE)/sv_inversion.bedpe.gz sv_inversion.bedpe.gz
	ln -s $(CURDATE)/sv_translocation.bedpe.gz sv_translocation.bedpe.gz
	
	ln -s $(CURDATE)/sv_deletion.bedpe.gz.tbi sv_deletion.bedpe.gz.tbi
	ln -s $(CURDATE)/sv_duplication.bedpe.gz.tbi sv_duplication.bedpe.gz.tbi
	ln -s $(CURDATE)/sv_insertion.bedpe.gz.tbi sv_insertion.bedpe.gz.tbi
	ln -s $(CURDATE)/sv_inversion.bedpe.gz.tbi sv_inversion.bedpe.gz.tbi
	ln -s $(CURDATE)/sv_translocation.bedpe.gz.tbi sv_translocation.bedpe.gz.tbi
	
	ln -s $(CURDATE)/sv_breakpoint_density.igv sv_breakpoint_density.igv
	ln -s $(CURDATE)/sv_breakpoint_density_*.igv .

export_somatic:
	mkdir -p $(CURDATE)
	$(NGSBITS)/NGSDExportAnnotationData -somatic $(CURDATE)/NGSD_somatic_unsorted.vcf > $(CURDATE)/NGSD_somatic.log 2>&1
	$(NGSBITS)/VcfStreamSort -in $(CURDATE)/NGSD_somatic_unsorted.vcf -out $(CURDATE)/NGSD_somatic.vcf
	bgzip -c $(CURDATE)/NGSD_somatic.vcf > $(CURDATE)/NGSD_somatic.vcf.gz
	tabix -p vcf $(CURDATE)/NGSD_somatic.vcf.gz
	#cleanup
	rm $(CURDATE)/NGSD_somatic_unsorted.vcf
	rm $(CURDATE)/NGSD_somatic.vcf
	#linking
	rm -f NGSD_somatic.vcf.gz
	rm -f NGSD_somatic.vcf.gz.tbi
	ln -s $(CURDATE)/NGSD_somatic.vcf.gz NGSD_somatic.vcf.gz
	ln -s $(CURDATE)/NGSD_somatic.vcf.gz.tbi NGSD_somatic.vcf.gz.tbi
