# preprocessing to viral genome db to remove redundancies and preferentially select our HPV-labeled viruses among sequence clusters.

cd-hit-est -i virus_db.fasta -o virus_db.fasta.cdhit -c .98 -d 0 

select_cluster_rep_HPV_pref.py > accs.selected

acc_to_FASTA_file.pl accs.selected virus_db.fasta > virus_db.nr.fasta



#---------------------------------------------------------------------
# Background Info: How the virus database fasta was originally compiled

(as posted here: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/wiki/Human-Virus-Database-Prep#background-info--how-the-virus-database-fasta-was-originally-compiled )
    
Nothing to do here wrt setting up CTAT-VIF... just info on how the above virus fasta file was constructed.

Below is documentation on how we created the virus database that's used with CTAT-VirusIntegrationFinder (CTAT-VIF).

Human viruses were downloaded from http://www.virusite.org

The list of human viruses were downloaded with parameter setting 'group=human' as file 'human_viruses.list.csv'.

All virus sequences from virusite.org were downloaded as 'genomes.fasta'

and then the subset of human viruses were extracted via:

CTAT-VirusIntegrationFinder/util/misc/extract_human_viruses.py | \
    perl -lane 's/_complete_(sequence|genome)//; s/refseq_//; print;'\
    >  human_viruses.fasta 
We prepended 143 HPV sequences to this fasta file as obtained from collaborating researchers.

To exclude additional occurrences of non-unique virus entries, we removed the redundant entries using cd-hit, and preferentially retained our HPV-labeled entries over the virussite entries.

cd-hit-est -i virus_db.fasta -o virus_db.fasta.cdhit -c .98 -d 0 

select_cluster_rep_HPV_pref.py > accs.selected

acc_to_FASTA_file.pl accs.selected virus_db.fasta > virus_db.nr.fasta
This final virus fasta file was named as virus_db.nr.fasta and is included as part of supplementary data resources to be used with CTAT-VIF.
