
CWD = $(CURDIR)


test:
	../ctat-vif \
		--left reads_1.fastq.gz \
		--right reads_2.fastq.gz \
		--sample_id test \
        --no_remove_duplicates \
        --genome_lib_dir ${CTAT_GENOME_LIB} \
		-O VIF.full.no_remove_dups.outdir

all: test star_init_only remove_dups novirus test_SE redherrings




test_kickstart_p2:

	../ctat-vif \
		--kickstart_p2_hg_unmapped_left_fq $(CWD)/__testdata/test.hgOnly.Unmapped.out.mate1.gz \
		--kickstart_p2_hg_unmapped_right_fq $(CWD)/__testdata/test.hgOnly.Unmapped.out.mate2.gz \
		--sample_id test \
        --no_remove_duplicates \
        --genome_lib_dir ${CTAT_GENOME_LIB} \
		-O VIF.full.no_remove_dups.outdir.kickstart_p2


test_kickstart_p3:

	../ctat-vif \
		--kickstart_p2_hg_unmapped_left_fq $(CWD)/__testdata/test.hgOnly.Unmapped.out.mate1.gz \
		--kickstart_p2_hg_unmapped_right_fq $(CWD)/__testdata/test.hgOnly.Unmapped.out.mate2.gz \
        --kickstart_p3_human_virus_chimJ $(CWD)/__testdata/human_virus_chimJ.tsv.gz \
        --kickstart_p3_human_virus_bam $(CWD)/__testdata/test.hgPlusVirus.Aligned.sortedByCoord.out.bam \
        --kickstart_p3_human_virus_bai $(CWD)/__testdata/test.hgPlusVirus.Aligned.sortedByCoord.out.bam.bai \
		--sample_id test \
        --no_remove_duplicates \
        --genome_lib_dir ${CTAT_GENOME_LIB} \
		-O VIF.full.no_remove_dups.outdir.kickstart_p3



test_SE:
	../ctat-vif \
		--left reads_1.fastq.gz \
		--sample_id test \
        --no_remove_duplicates \
        --genome_lib_dir ${CTAT_GENOME_LIB} \
		-O VIF.full.no_remove_dups.SE.outdir


star_init_only:
	../ctat-vif \
		--left reads_1.fastq.gz \
		--right reads_2.fastq.gz \
		--sample_id test \
		-O VIF.star_init_only.outdir \
		--star_init_only

remove_dups:
	../ctat-vif \
		--left reads_1.fastq.gz \
		--right reads_2.fastq.gz \
		--sample_id test \
		-O VIF.full.remove_dups.outdir


docker:
	docker run --rm -it \
		-v `pwd`:/data \
		-v ${CTAT_GENOME_LIB}:/ctat_genome_lib trinityctat/ctat_vif:latest \
		/usr/local/src/CTAT-VirusIntegrationFinder/ctat-vif \
		--left /data/reads_1.fastq.gz \
		--right /data/reads_2.fastq.gz \
		--sample_id test \
		--genome_lib_dir /ctat_genome_lib \
		-O /data/VIF.docker.outdir \
		--cpu 10


singularity:
	singularity exec -e -H `pwd` -B ${CTAT_GENOME_LIB} ../Docker/ctat_vif.*.simg ctat-vif \
		--left reads_1.fastq.gz \
		--right reads_2.fastq.gz \
		--sample_id test \
		--genome_lib_dir ${CTAT_GENOME_LIB}  \
		-O VIF.singularity.outdir \
		--cpu 10



redherrings:
	../ctat-vif \
		--left red_herrings_1.fastq.gz \
		--right red_herrings_2.fastq.gz \
		--sample_id redherr \
        --genome_lib_dir ${CTAT_GENOME_LIB} \
		-O VIF.redherrings.outdir





novirus:
	../ctat-vif \
		--left test.novirus.reads_1.fastq.gz \
		--right test.novirus.reads_2.fastq.gz \
		--sample_id test \
		-O VIF.novirus.outdir




LRsplice:
	../util/ctat-VIF.longrange_virus_splice.py \
		--left reads_1.fastq.gz \
		--right reads_2.fastq.gz \
		--virus_insertions_tsv __testdata/test.GRCH38.insertion_candidates.tsv \
		--sample_id test \
		--viral_fasta  ${CTAT_GENOME_LIB}/viruses/HPVs_db.fasta --viral_db_gtf  ${CTAT_GENOME_LIB}/viruses/HPVs_db.annots.gtf --flank 1 



clean:
	rm -rf ./VIF.*.outdir


