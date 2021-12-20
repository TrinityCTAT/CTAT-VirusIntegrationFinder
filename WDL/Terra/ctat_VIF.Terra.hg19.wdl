version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.0.1/WDL/ctat_VIF.wdl" as ctat_VIF_wf


workflow ctat_VIF_Terra_hg19 {

  input {
    String sample_id
    File left
    File? right

    }
  
  call ctat_VIF_wf {
    input:     
      sample_id = sample_id,
    left = left,
    right = right,
    viral_fasta = "




   }

}
