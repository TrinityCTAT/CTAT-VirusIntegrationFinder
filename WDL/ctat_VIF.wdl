version 1.0

workflow ctat_vif {
    input {
        File left
        File? right

        File fasta
        File gtf
        File viral_fasta
        File? viral_gtf

        Boolean star_init_only = false
        Boolean remove_duplicates = false
        Int min_reads = 0

        File? star_reference
        String? star_reference_dir

        String util_dir = "/usr/local/bin/CTAT-VirusIntegrationFinder/util"
        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        Boolean star_use_ssd = false
        Int star_cpu = 12
        Float star_memory = 75
        Int preemptible = 2
        String docker = "trinityctat/ctat_vif:0.1.0"
    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}

        min_reads:{help:"Filter insertion sites candidates that do not have at least 'min_reads'"}
        remove_duplicates:{help:"Remove duplicate alignments"}

        viral_fasta:{help:"Viral fasta"}

        star_reference:{help:"STAR index archive"}
        star_reference_dir:{help:"STAR directory (for non-Terra use)"}
        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        util_dir:{help:"Path to util directory"}


        star_init_only:{help:"Only perform initial STAR chimeric junction analysis"}
        docker:{help:"Docker image"}
    }

    output {
        File star_bam = STAR.bam
        File star_bam_index = STAR.bai
        File star_output_log_final = STAR.output_log_final
        File star_output_SJ = STAR.output_SJ
        File? star_chimeric_junction = STAR.chimeric_junction

        File? remove_duplicates_bam = RemoveDuplicates.bam
        File? remove_duplicates_bam_index = RemoveDuplicates.bai

        File insertion_site_candidates_full = InsertionSiteCandidates.full
        File insertion_site_candidates_abridged = InsertionSiteCandidates.abridged
        File? insertion_site_candidates_abridged_filtered = InsertionSiteCandidates.abridged_filtered
        File insertion_site_candidates_abridged_detailed = InsertionSiteCandidates.abridged_detailed

        File genome_abundance_plot = GenomeAbundancePlot.plot

        File virus_coverage_read_counts_summary = TopVirusCoverage.read_counts_summary
        File virus_coverage_read_counts_image = TopVirusCoverage.read_counts_image
        File virus_coverage_read_counts_log_image = TopVirusCoverage.read_counts_log_image
        Array[File] virus_coverage_virus_images = TopVirusCoverage.virus_images

        File igv_virus_report_html = IGVVirusReport.html

        File? fasta_extract = ExtractChimericGenomicTargets.fasta_extract
        File? gtf_extract = ExtractChimericGenomicTargets.gtf_extract

        File? star2_bam = STAR2.bam
        File? star2_bam_index = STAR2.bai
        File? star2_output_log_final = STAR2.output_log_final
        File? star2_output_SJ = STAR2.output_SJ

        File? remove_duplicates2_bam = RemoveDuplicates2.bam
        File? remove_duplicates2_bam_index = RemoveDuplicates2.bai

        File? evidence_counts =  ChimericContigEvidenceAnalyzer.evidence_counts
        File? evidence_bam = ChimericContigEvidenceAnalyzer.evidence_bam
        File? evidence_bai = ChimericContigEvidenceAnalyzer.evidence_bai

        File? refined_counts = RefineVIFOutput.refined_counts

        File? genome_abundance_refined_plot = GenomeAbundancePlot2.plot

        File? igv_report_html = IGVReport.html
    }
    call STAR {
        input:
            util_dir=util_dir,
            fastq1=left,
            fastq2=right,
            base_name="out",
            star_reference=star_reference,
            star_reference_dir=star_reference_dir,
            viral_fasta=viral_fasta,
            viral_gtf=viral_gtf,
            disable_chimeras=false,
            extra_disk_space = star_extra_disk_space,
            disk_space_multiplier = star_fastq_disk_space_multiplier,
            memory = star_memory,
            use_ssd = star_use_ssd,
            cpu = star_cpu,
            docker = docker,
            preemptible = preemptible
    }
    if(remove_duplicates) {
        call RemoveDuplicates {
            input:
                input_bam=STAR.bam,
                input_bai=STAR.bai,
                output_bam = "Aligned.sortedByCoord.out.rm.dups.bam",
                util_dir=util_dir,
                cpu=1,
                preemptible=preemptible,
                memory="2G",
                docker=docker,
                extra_disk_space=1,
                disk_space_multiplier=2
        }
    }
    call InsertionSiteCandidates {
        input:
            chimeric_junction=select_first([STAR.chimeric_junction]),
            viral_fasta=viral_fasta,
            remove_duplicates=remove_duplicates,
            util_dir=util_dir,
            min_reads=min_reads,
            cpu=1,
            preemptible=preemptible,
            memory="2G",
            docker=docker,
    }
    File insertion_site_candidates_output = select_first([InsertionSiteCandidates.abridged_filtered, InsertionSiteCandidates.abridged])

    call GenomeAbundancePlot {
        input:
            counts=insertion_site_candidates_output,
            output_name="vif.prelim.genome_plot",
            title="Preliminary Genome Wide Abundance",
            util_dir=util_dir,
            preemptible=preemptible,
            docker=docker
    }
    File aligned_bam = select_first([RemoveDuplicates.bam, STAR.bam])
    File aligned_bai = select_first([RemoveDuplicates.bai, STAR.bai])
    call TopVirusCoverage {
        input:
            chimeric_events=insertion_site_candidates_output,
            bam=aligned_bam,
            bai=aligned_bai,
            util_dir=util_dir,
            preemptible=preemptible,
            docker=docker
    }
    call IGVVirusReport {
        input:
            bam=aligned_bam,
            bai=aligned_bai,
            viral_fasta=viral_fasta,
            read_counts_summary=TopVirusCoverage.read_counts_summary,
            util_dir=util_dir,
            preemptible=preemptible,
            docker=docker
    }

    if(!star_init_only) {
        call ExtractChimericGenomicTargets {
            input:
                bam=aligned_bam,
                bai=aligned_bai,
                fasta=fasta,
                viral_fasta=viral_fasta,
                insertion_site_candidates_abridged=insertion_site_candidates_output,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
        call STAR as STAR2 {
            input:
                util_dir=util_dir,
                fastq1=left,
                fastq2=right,
                base_name="out2",
                star_reference=star_reference,
                star_reference_dir=star_reference_dir,
                viral_fasta=ExtractChimericGenomicTargets.fasta_extract,
                disable_chimeras=true,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                docker = docker,
                preemptible = preemptible
        }
        if(remove_duplicates) {
            call RemoveDuplicates as RemoveDuplicates2 {
                input:
                    input_bam=STAR2.bam,
                    input_bai=STAR2.bai,
                    output_bam = "Aligned2.sortedByCoord.out.rm.dups.bam",
                    util_dir=util_dir,
                    cpu=1,
                    preemptible=preemptible,
                    memory="2G",
                    docker=docker,
                    extra_disk_space=1,
                    disk_space_multiplier=2
            }
        }
        File aligned_bam2 = select_first([RemoveDuplicates2.bam, STAR2.bam])
        File aligned_bai2 = select_first([RemoveDuplicates2.bai, STAR2.bai])

        call ChimericContigEvidenceAnalyzer {
            input:
                bam=aligned_bam2,
                bai=aligned_bai2,
                gtf=ExtractChimericGenomicTargets.gtf_extract,
                min_reads=min_reads,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
        call RefineVIFOutput {
            input:
                prelim_counts=insertion_site_candidates_output,
                vif_counts=ChimericContigEvidenceAnalyzer.evidence_counts,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
        call GenomeAbundancePlot as GenomeAbundancePlot2{
            input:
                counts=RefineVIFOutput.refined_counts,
                output_name="vif.genome_plot",
                title="Genome Wide Abundance",
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
        call IGVReport {
            input:
                summary_results_tsv=RefineVIFOutput.refined_counts,
                alignment_bam=ChimericContigEvidenceAnalyzer.evidence_bam,
                alignment_bai=ChimericContigEvidenceAnalyzer.evidence_bai,
                chim_targets_gtf=ExtractChimericGenomicTargets.gtf_extract,
                chim_targets_fasta=ExtractChimericGenomicTargets.fasta_extract,
                gtf=gtf,
                images=[GenomeAbundancePlot.plot, GenomeAbundancePlot2.plot, TopVirusCoverage.read_counts_image, TopVirusCoverage.read_counts_log_image],
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
    }
}


task STAR {
    input {
        String util_dir
        File fastq1
        File? fastq2
        File? star_reference
        String? star_reference_dir
        File viral_fasta
        File? viral_gtf
        Boolean disable_chimeras
        Int max_mate_dist = 100000

        Float extra_disk_space
        Float disk_space_multiplier
        Boolean use_ssd
        Int cpu
        Int preemptible
        String memory
        String docker
        String base_name

    }

    Boolean is_gzip = sub(select_first([fastq1]), "^.+\\.(gz)$", "GZ") == "GZ"
    command <<<

        genomeDir="~{star_reference}"
        if [ "$genomeDir" == "" ]; then
            genomeDir="~{star_reference_dir}"
        fi

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            tar xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN ~{cpu} \
        --readFilesIn ~{fastq1} ~{fastq2} \
        ~{true='--readFilesCommand "gunzip -c"' false='' is_gzip} \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ~{base_name}. \
        --twopassMode Basic \
        --alignSJDBoverhangMin 10 \
        --genomeSuffixLengthMax 10000 \
        --limitBAMsortRAM 47271261705 \
        --alignInsertionFlush Right \
        --alignMatesGapMax ~{max_mate_dist} \
        --alignIntronMax  ~{max_mate_dist} \
        --genomeFastaFiles ~{viral_fasta} \
        --outSAMfilter KeepAllAddedReferences \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --scoreGapNoncan -6 \
        ~{true='' false='--chimJunctionOverhangMin 12 --chimSegmentMin 12 --chimSegmentReadGapMax 3' disable_chimeras} \
        ~{"--sjdbOverhang 150 --sjdbGTFfile " + viral_gtf}

        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"
    >>>

    output {
        File bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File bai = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"
        File output_log_final = "~{base_name}.Log.final.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        File? chimeric_junction = "~{base_name}.Chimeric.out.junction"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*disk_space_multiplier + size(fastq2, "GB") * disk_space_multiplier + size(star_reference, "GB")*8 + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }
}

task RemoveDuplicates {
    input {
        File input_bam
        File input_bai
        String util_dir
        Int cpu
        Int preemptible
        String memory
        String docker
        Float extra_disk_space
        Float disk_space_multiplier
        String output_bam
    }
    command <<<
        ~{util_dir}/bam_mark_duplicates.py \
        -i ~{input_bam} \
        -o ~{output_bam} \
        -r

        samtools index ~{output_bam}
    >>>

    output {
        File bam = "~{output_bam}"
        File bai = "~{output_bam}.bai"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(input_bam, "GB")*disk_space_multiplier + extra_disk_space) + " HDD"
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }
}

task InsertionSiteCandidates {
    input {
        File chimeric_junction
        File viral_fasta
        Boolean remove_duplicates
        String util_dir
        Int cpu
        Int preemptible
        String memory
        String docker
        Int min_reads

    }
    String prefix = "vif.prelim"

    command <<<
        ~{util_dir}/chimJ_to_virus_insertion_candidate_sites.py \
        --chimJ ~{chimeric_junction} \
        --patch_db_fasta ~{viral_fasta} \
        --output_prefix ~{prefix} \
        ~{true='--remove_duplicates' false='' remove_duplicates}

        python <<CODE
        import pandas as pd
        min_reads = ~{min_reads}
        if min_reads > 0:
            df = pd.read_csv("~{prefix}.abridged.tsv", sep='\t')
            df = df[df['total'] >= min_reads]
            df.to_csv("~{prefix}.abridged.filtered.tsv", sep='\t', index=False)
        CODE
    >>>

    output {
        File full = "~{prefix}.full.tsv"
        File abridged = "~{prefix}.abridged.tsv"
        File? abridged_filtered = "~{prefix}.abridged.filtered.tsv"
        File abridged_detailed = "~{prefix}.abridged.detailed.tsv"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(viral_fasta, "GB") + size(chimeric_junction, "GB")*3) + " HDD"
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }
}

task TopVirusCoverage {
    input {
        File chimeric_events
        String util_dir
        File bam
        File bai
        Int preemptible
        String docker
    }
    String prefix = "vif"

    command <<<
        ~{util_dir}/plot_top_virus_coverage.Rscript \
        --vif_report ~{chimeric_events} \
        --bam ~{bam} \
        --output_prefix ~{prefix}
    >>>

    output {
        File read_counts_summary = "~{prefix}.virus_read_counts_summary.tsv"
        File read_counts_image = "~{prefix}.virus_read_counts.png"
        File read_counts_log_image = "~{prefix}.virus_read_counts_log.png"
        Array[File] virus_images = glob("~{prefix}.virus_coverage_*.png")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(chimeric_events, "GB") + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}



task ExtractChimericGenomicTargets {
    input {
        File bam
        File bai
        File fasta
        File viral_fasta
        File insertion_site_candidates_abridged
        String util_dir
        Int preemptible
        String docker
    }

    command <<<

        ~{util_dir}/extract_chimeric_genomic_targets.py \
        --fasta ~{fasta} \
        --patch_db_fasta ~{viral_fasta} \
        --output_prefix vif.extract \
        --chim_events ~{insertion_site_candidates_abridged} \
        --pad_region_length 1000
    >>>

    output {
        File fasta_extract = "vif.extract.fasta"
        File gtf_extract = "vif.extract.gtf"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(viral_fasta, "GB")*2 + size(fasta, "GB")*2 + size(bam, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task ChimericContigEvidenceAnalyzer {
    input {
        File bam
        File bai
        File gtf

        String util_dir
        Int preemptible
        String docker
        Int min_reads
    }

    command <<<

        ~{util_dir}/chimeric_contig_evidence_analyzer.py \
        --patch_db_bam ~{bam} \
        --patch_db_gtf ~{gtf} \
        --output_prefix vif

        samtools index vif.evidence.bam
    >>>

    output {
        File evidence_counts =  "vif.evidence_counts.tsv"
        File evidence_bam ="vif.evidence.bam"
        File evidence_bai ="vif.evidence.bam.bai"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(bam, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task RefineVIFOutput {
    input {
        File prelim_counts # InsertionSiteCandidates.abridged
        File vif_counts # ChimericContigEvidenceAnalyzer.evidence_counts

        String util_dir
        Int preemptible
        String docker
    }

    command <<<

        ~{util_dir}/refine_VIF_output.Rscript \
        --prelim_counts ~{prelim_counts} \
        --vif_counts ~{vif_counts} \
        --output vif.refined.tsv
    >>>

    output {
        File refined_counts =  "vif.refined.tsv"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(prelim_counts, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task GenomeAbundancePlot {
    input {
        File counts
        String output_name
        String title
        String util_dir
        Int preemptible
        String docker
    }

    command <<<
        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
        --vif_report ~{counts} \
        --title "~{title}" \
        --output_png ~{output_name}.png
    >>>

    output {
        File plot =  "~{output_name}.png"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(counts, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task IGVVirusReport {
    input {
        File bam
        File bai
        File viral_fasta
        File read_counts_summary
        String util_dir
        Int preemptible
        String docker
    }
    Int max_coverage = 100

    command <<<

        ~{util_dir}/create_insertion_site_inspector_js.py \
        --VIF_summary_tsv ~{read_counts_summary} \
        --json_outfile vif.virus.json

        # make bed for igvjs
        ~{util_dir}/create_virus_bed.py ~{read_counts_summary} vif.virus.bed

        # prep for making the report
        ~{util_dir}/bamsifter/bamsifter \
        -c ~{max_coverage} \
        -o vif.virus.reads.bam \
        ~{bam}

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        ln -sf ~{viral_fasta} vif.virus.fa

        # generate the html
        ~{util_dir}/make_VIF_igvjs_html.py \
        --html_template ~{util_dir}/resources/igvjs_VIF.html \
        --fusions_json vif.virus.json \
        --input_file_prefix vif.virus \
        --html_output vif.virus.html
    >>>

    output {
        File html = "vif.virus.html"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(bam, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task IGVReport {
    input {
        File summary_results_tsv
        File alignment_bam
        File alignment_bai
        File chim_targets_gtf
        File chim_targets_fasta
        File gtf
        Array[File] images
        String util_dir
        Int preemptible
        String docker
    }
    Int max_coverage = 100

    command <<<
        set -ex

        ~{util_dir}/find_closest.py \
        -i ~{summary_results_tsv} \
        -o summary_results_tsv_with_genes.tsv \
        --gtf ~{gtf}

        ~{util_dir}/create_insertion_site_inspector_js.py \
        --VIF_summary_tsv summary_results_tsv_with_genes.tsv \
        --json_outfile igv.json

        # make bed for igvjs
        ~{util_dir}/region_gtf_to_bed.py ~{chim_targets_gtf} > vif.bed

        # prep for making the report
        ~{util_dir}/bamsifter/bamsifter \
        -c ~{max_coverage} \
        -o vif.reads.bam \
        ~{alignment_bam}

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        ln -sf ~{chim_targets_fasta} vif.fa

        ~{util_dir}/make_VIF_igvjs_html.py \
        --html_template ~{util_dir}/resources/igvjs_VIF.html \
        --fusions_json igv.json \
        --input_file_prefix vif \
        --html_output igv.tmp.html

        # generate the final report
        ~{util_dir}/add_to_html.py \
        --html igv.tmp.html \
        --out igv.html \
        --image ~{sep=' --image ' images}
    >>>

    output {
        File html = "igv.html"
    }
    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(alignment_bam, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}
