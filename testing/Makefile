


test:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.full.no_remove_dups.outdir

star_init_only:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.star_init_only.outdir --star_init_only

remove_dups:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.full.remove_dups.outdir --remove_duplicates



singularity:
	singularity exec -e -H `pwd` -B {$CTAT_GENOME_LIB} ../Docker/ctat_vif.*.simg ctat-VIF.py --viral_db_fasta HPV16.fa --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --genome_lib_dir ${CTAT_GENOME_LIB}  -O VIF.singularity.outdir --CPU 10



clean:
	rm -rf ./VIF.*.outdir


