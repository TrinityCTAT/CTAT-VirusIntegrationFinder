# preprocessing to viral genome db to remove redundancies and preferentially select our HPV-labeled viruses among sequence clusters.

cd-hit-est -i virus_db.fasta -o virus_db.fasta.cdhit -c .98 -d 0 

select_cluster_rep_HPV_pref.py > accs.selected

acc_to_FASTA_file.pl accs.selected virus_db.fasta > virus_db.nr.fasta

