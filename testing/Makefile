


test:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.full.no_remove_dups.outdir

all: test star_init_only remove_dups novirus


star_init_only:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.star_init_only.outdir --star_init_only

remove_dups:
	../ctat-VIF.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.full.remove_dups.outdir --remove_duplicates



singularity:
	singularity exec -e -H `pwd` -B {$CTAT_GENOME_LIB} ../Docker/ctat_vif.*.simg ctat-VIF.py --viral_db_fasta HPV16.fa --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --genome_lib_dir ${CTAT_GENOME_LIB}  -O VIF.singularity.outdir --CPU 10


novirus:
	../ctat-VIF.py --left_fq test.novirus.reads_1.fastq.gz --right_fq test.novirus.reads_2.fastq.gz --viral_db_fasta HPV16.fa -O VIF.novirus.outdir




LRsplice:
	../util/ctat-VIF.longrange_virus_splice.py --left_fq V42_hpv16chim_1.fastq.gz --right_fq V42_hpv16chim_2.fastq.gz --virus_insertions_tsv __testdata/test.GRCH38.insertion_candidates.tsv --viral_db_fasta  /seq/RNASEQ/HPV_EinsteinMed/HPV_Genomes/All_HPV_Types/HPVs_db.fasta --viral_db_gtf  /seq/RNASEQ/HPV_EinsteinMed/HPV_Genomes/All_HPV_Types/HPVs_db.annots.gtf --flank 1 



clean:
	rm -rf ./VIF.*.outdir


