version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.0.1c/WDL/Terra/ctat_VIF.Terra.wdl" as ctat_VIF_Terra
#import "ctat_VIF.Terra.wdl" as ctat_VIF_Terra


workflow ctat_VIF_Terra_hg38 {

  input {
    String sample_id
    File? left
    File? right
    File? drs_path_fastqs
    String docker = "trinityctat/ctat_vif:1.0.1"
    Int preemptible = 0

    CTAT_VIF_config pipe_inputs_config = {
      "ref_genome_fasta" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa",
      "ref_genome_gtf" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_annot.gtf",
      "viral_fasta" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/Virus_db_Dec072021.nonUnq100bMsk.filt90.fasta",
      "star_index_human_only" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa.star.idx.tar",
      "star_index_human_plus_virus" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/GRCh38_plus_viraldb_Dec092021.fasta.star.idx.tar",
      "NULL_file" : "gs://ctat_genome_libs/null" 
    }

  }
  
  call ctat_VIF_Terra.ctat_VIF_Terra {
    input:     
      sample_id = sample_id,
      left = left,
      right = right,
      drs_path_fastqs = drs_path_fastqs,
      docker = docker,
      pipe_inputs_config = pipe_inputs_config,
      preemptible = preemptible
   }

}

