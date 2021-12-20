version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.0.1/WDL/ctat_VIF.wdl" as ctat_VIF_wf


struct CTAT_VIF_config {

  File ref_genome_fasta
  File ref_genome_gtf
  File viral_fasta
  File star_index_human_only
  File star_index_human_plus_virus

}


workflow ctat_VIF_Terra {

  input {
    String sample_id
    File left
    File? right
    String docker = "trinityctat/ctat_vif:1.0.1"
    CTAT_VIF_config pipe_inputs_config

    }
  
  call ctat_VIF_wf.ctat_vif {
    input:     
      sample_id = sample_id,
      left = left,
      right = right,
      docker = docker,
    
      ref_genome_fasta = pipe_inputs_config.ref_genome_fasta,
      ref_genome_gtf = pipe_inputs_config.ref_genome_gtf,
      viral_fasta = pipe_inputs_config.viral_fasta,
      star_index_human_only = pipe_inputs_config.star_index_human_only,
      star_index_human_plus_virus = pipe_inputs_config.star_index_human_plus_virus

   }

}

