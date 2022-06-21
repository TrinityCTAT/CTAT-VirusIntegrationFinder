version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.3.1/WDL/Terra/ctat_VIF.Terra.wdl" as ctat_VIF_Terra
#import "ctat_VIF.Terra.wdl" as ctat_VIF_Terra


workflow ctat_VIF_Terra_hg19 {

  input {
    String sample_id
    File left
    File? right
    File? drs_path_fastqs
    Boolean clean_reads = true
    Int max_hits = 50
    String docker = "trinityctat/ctat_vif:latest"
    Int preemptible = 0
    
    CTAT_VIF_config pipe_inputs_config = {
      "ref_genome_fasta" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.fa",
      "ref_genome_gtf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_annot.gtf",
      "viral_fasta" : "gs://ctat_genome_libs/GRCh38_gencode_v22/06-06-2022/virus_db.fasta",
      "star_index_human_only" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.fa.star.idx.tar",
      "star_index_human_plus_virus" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/06-13-2022/hg_plus_viraldb.fasta.star.idx.tar",
      "NULL_file" : "gs://ctat_genome_libs/null"
    }
    
  }
  
  call ctat_VIF_Terra.ctat_VIF_Terra {
    input:     
      sample_id = sample_id,
      left = left,
      right = right,
      clean_reads = clean_reads,
      max_hits = max_hits,
      drs_path_fastqs = drs_path_fastqs,
      docker = docker,
      pipe_inputs_config = pipe_inputs_config,
      preemptible = preemptible
   }

}

