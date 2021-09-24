version 1.0


workflow ctat_vif {
    input {
        File? left
        File? right

        File fasta
        File gtf

        File viral_fasta
        File? viral_gtf

        Boolean star_init_only = false
        Boolean remove_duplicates = false
        Boolean generate_reports = true
        Int min_reads = 0

        File? star_reference
        String? star_reference_dir

        String igv_virus_reports_memory = "14GB"
        String igv_reports_memory = "1GB"

        String util_dir = "/usr/local/src/CTAT-VirusIntegrationFinder/util"
        String picard = "/usr/local/src/picard.jar"
        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        String star_index_memory = "50G"
        Int sjdb_overhang = 150

        String? sample_id

        Boolean autodetect_cpu = true # auto-detect number of cpus for STAR as # of requested CPUs might not equal actual CPUs, depending on memory
        Boolean star_use_ssd = false
        Int star_cpu = 12
        Float star_init_memory = 45
        String star_init_two_pass_mode = "Basic"

        Float star_memory = 45
        String star_two_pass_mode = "Basic" # or None
        Int preemptible = 2
        String docker = "trinityctat/ctat_vif:0.1.0"

        # run stage 2 only
        File? bam
        File? bam_index
        File? insertion_site_candidates
    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        min_reads:{help:"Filter insertion sites candidates that do not have at least 'min_reads'"}
        remove_duplicates:{help:"Remove duplicate alignments"}
        fasta:{help:"Host fasta"}
        gtf:{help:"Host annotations GTF"}
        viral_fasta:{help:"Viral fasta"}
        bam:{help:"Previously aligned bam file"}
        bam_index:{help:"BAM index corresponding to bam file"}
        insertion_site_candidates:{help:"Previously generated candidates"}
        star_reference:{help:"STAR index archive containing both host and viral genomes"}
        star_reference_dir:{help:"STAR directory containing both host and viral genomes (for non-Terra use)"}
        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        util_dir:{help:"Path to util directory (for non-Docker use)"}
        star_init_only:{help:"Only perform initial STAR chimeric junction analysis"}
        docker:{help:"Docker image"}
    }

    output {
        File? star_bam = STAR.bam
        File? star_bam_index = STAR.bai
        File? star_output_log_final = STAR.output_log_final
        File? star_output_SJ = STAR.output_SJ
        File? star_chimeric_junction = STAR.chimeric_junction

        File? remove_duplicates_bam = RemoveDuplicates.bam
        File? remove_duplicates_bam_index = RemoveDuplicates.bai

        File? insertion_site_candidates_full = InsertionSiteCandidates.full
        File? insertion_site_candidates_abridged = InsertionSiteCandidates.abridged
        File? insertion_site_candidates_abridged_filtered = InsertionSiteCandidates.abridged_filtered
        File? insertion_site_candidates_abridged_detailed = InsertionSiteCandidates.abridged_detailed

        File? genome_abundance_plot = VirusReport.genome_abundance_plot
        File? virus_coverage_read_counts_summary = VirusReport.read_counts_summary
        File? virus_coverage_read_counts_image = VirusReport.read_counts_image
        File? virus_coverage_read_counts_log_image = VirusReport.read_counts_log_image
        Array[File]? virus_coverage_virus_images = VirusReport.virus_images
        File? igv_virus_report_html = VirusReport.html

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

        File? refined_counts = SummaryReport.refined_counts
        File? genome_abundance_refined_plot = SummaryReport.genome_abundance_plot
        File? igv_report_html = SummaryReport.html

    }

    Boolean create_star_index = !defined(star_reference) && !defined(star_reference_dir)
    if(create_star_index) {
        call STARIndex {
            input:
                fasta=fasta,
                gtf=gtf,
                sjdb_overhang=sjdb_overhang,
                viral_fasta=viral_fasta,
                memory = star_index_memory,
                use_ssd = star_use_ssd,
                autodetect_cpu = autodetect_cpu,
                cpu = star_cpu,
                docker = docker,
                preemptible = preemptible
        }
    }

    File? star_reference_use = (if(create_star_index) then STARIndex.genome else star_reference)
    if(!defined(bam) && defined(left)) {
        call STAR {
            input:
                util_dir=util_dir,
                fastq1=select_first([left]),
                fastq2=right,
                two_pass_mode = star_init_two_pass_mode,
                base_name="out",
                star_reference=star_reference_use,
                star_reference_dir=star_reference_dir,
                disable_chimeras=false,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_init_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                autodetect_cpu = autodetect_cpu,
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
                preemptible=preemptible,
                docker=docker,
        }
        File insertion_site_candidates_output =  (if min_reads>0 then select_first([InsertionSiteCandidates.abridged_filtered]) else InsertionSiteCandidates.abridged)

        if(generate_reports) {
            call VirusReport {
                input:
                    bam=select_first([RemoveDuplicates.bam, STAR.bam]),
                    bai=select_first([RemoveDuplicates.bai, STAR.bai]),
                    viral_fasta=viral_fasta,
                    insertion_site_candidates=insertion_site_candidates_output,
                    util_dir=util_dir,
                    memory=igv_virus_reports_memory,
                    preemptible=preemptible,
                    docker=docker
            }
        }

    }

    if(!star_init_only) {
        File aligned_bam_use = select_first([bam, RemoveDuplicates.bam, STAR.bam])
        File aligned_bai_use= select_first([bam_index, RemoveDuplicates.bai, STAR.bai])
        File insertion_site_candidates_use = select_first([insertion_site_candidates, insertion_site_candidates_output])
        call ExtractChimericGenomicTargets {
            input:
                bam=aligned_bam_use,
                bai=aligned_bai_use,
                fasta=fasta,
                viral_fasta=viral_fasta,
                insertion_site_candidates_abridged=insertion_site_candidates_use,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker
        }
        call STAR as STAR2 {
            input:
                util_dir=util_dir,
                fastq1=select_first([left]),
                fastq2=right,
                two_pass_mode = star_two_pass_mode,
                base_name="out2",
                star_reference=star_reference_use,
                star_reference_dir=star_reference_dir,
                genome_fasta_file=ExtractChimericGenomicTargets.fasta_extract,
                disable_chimeras=true,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                autodetect_cpu = autodetect_cpu,
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

        if(generate_reports) {
            call SummaryReport {
                input:
                    prelim_counts=insertion_site_candidates_use,
                    vif_counts=ChimericContigEvidenceAnalyzer.evidence_counts,
                    alignment_bam=ChimericContigEvidenceAnalyzer.evidence_bam,
                    alignment_bai=ChimericContigEvidenceAnalyzer.evidence_bai,
                    chim_targets_gtf=ExtractChimericGenomicTargets.gtf_extract,
                    chim_targets_fasta=ExtractChimericGenomicTargets.fasta_extract,
                    gtf=gtf,
                    images=select_all([VirusReport.genome_abundance_plot, VirusReport.read_counts_image, VirusReport.read_counts_log_image]),
                    util_dir=util_dir,
                    memory=igv_reports_memory,
                    preemptible=preemptible,
                    docker=docker
            }
        }
    }
}

task STARIndex {

    input {
        File fasta
        File gtf
        File viral_fasta
        Int cpu
        Int preemptible
        String memory
        String docker
        Boolean use_ssd
        Int sjdb_overhang
        Boolean autodetect_cpu
    }

    Float extra_disk_space = 100
    Float disk_space_multiplier = 10

    command <<<
        set -e

        cpu=~{cpu}
        if [ "~{autodetect_cpu}" == "true" ]; then
            cpu=$(nproc)
        fi

        mkdir genome_dir

        STAR \
        --runThreadN $cpu \
        --runMode genomeGenerate \
        --genomeDir genome_dir \
        --genomeFastaFiles ~{fasta} ~{viral_fasta} \
        ~{"--sjdbOverhang " + sjdb_overhang + " --sjdbGTFfile " + gtf}

        tar -I pigz -cf STAR.tar.gz genome_dir
    >>>

    output {
        File genome = "STAR.tar.gz"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fasta, "GB")*disk_space_multiplier + size(viral_fasta, "GB") * disk_space_multiplier + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory
    }
}

task STAR {
    input {
        String util_dir
        File fastq1
        File? fastq2
        File? star_reference
        String? star_reference_dir
        Boolean disable_chimeras
        Float extra_disk_space
        Float disk_space_multiplier
        Boolean use_ssd
        Int cpu
        Int preemptible
        Float memory
        String docker
        String base_name
        String two_pass_mode
        Boolean autodetect_cpu
        File? genome_fasta_file
    }
    Int max_mate_dist = 100000
    Boolean has_genome_fasta_file = defined(genome_fasta_file)

    command <<<
        set -e

        cpu=~{cpu}
        genomeDir="~{star_reference}"
        fastqs="~{fastq1} ~{fastq2}"
        readFilesCommand=""
        if [[ "~{fastq1}" = *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
        fi
        if [ "~{autodetect_cpu}" == "true" ]; then
            cpu=$(nproc)
        fi

        if [ "$genomeDir" == "" ]; then
            genomeDir="~{star_reference_dir}"
        fi

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            compress="pigz"

            if [[ $genomeDir == *.bz2 ]] ; then
                compress="pbzip2"
            fi
            tar -I $compress -xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            readFilesCommand=""
            if [[ "$fastqs" = *.gz ]] ; then
                readFilesCommand="--readFilesCommand \"gunzip -c\""
            fi
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN $cpu \
        --readFilesIn $fastqs \
        $readFilesCommand \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ~{base_name}. \
        ~{"--twopassMode " + two_pass_mode} \
        --alignSJDBoverhangMin 10 \
        --genomeSuffixLengthMax 10000 \
        --limitBAMsortRAM 47271261705 \
        --alignInsertionFlush Right \
        ~{"--genomeFastaFiles " + genome_fasta_file} \
        ~{true='--outSAMfilter KeepAllAddedReferences' false='' has_genome_fasta_file} \
        --alignMatesGapMax ~{max_mate_dist} \
        --alignIntronMax  ~{max_mate_dist} \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --scoreGapNoncan -6 \
        ~{true='' false='--chimJunctionOverhangMin 12 --chimSegmentMin 12 --chimSegmentReadGapMax 3' disable_chimeras}

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
        set -e

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
        memory: memory
    }
}

task InsertionSiteCandidates {
    input {
        File chimeric_junction
        File viral_fasta
        Boolean remove_duplicates
        String util_dir
        Int preemptible
        String docker
        Int min_reads

    }
    String prefix = "vif.prelim"

    command <<<
        set -e

        ~{util_dir}/chimJ_to_virus_insertion_candidate_sites.py \
        --chimJ ~{chimeric_junction} \
        --patch_db_fasta ~{viral_fasta} \
        --output_prefix ~{prefix} \
        ~{true='--remove_duplicates' false='' remove_duplicates}

        python <<CODE
        min_reads = ~{min_reads}
        if min_reads > 0:
            import pandas as pd
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
        cpu: 1
        memory: "1GB"
    }
}

#task TopVirusCoverage {
#    input {
#        File chimeric_events
#        String util_dir
#        File bam
#        File bai
#        Int preemptible
#        String docker
#    }
#    String prefix = "vif"
#
#    command <<<
#        set -e
#
#        ~{util_dir}/plot_top_virus_coverage.Rscript \
#        --vif_report ~{chimeric_events} \
#        --bam ~{bam} \
#        --output_prefix ~{prefix}
#    >>>
#
#    output {
#        File read_counts_summary = "~{prefix}.virus_read_counts_summary.tsv"
#        File read_counts_image = "~{prefix}.virus_read_counts.png"
#        File read_counts_log_image = "~{prefix}.virus_read_counts_log.png"
#        Array[File] virus_images = glob("~{prefix}.virus_coverage_*.png")
#    }
#
#    runtime {
#        preemptible: preemptible
#        disks: "local-disk " + ceil(size(bam, "GB") + size(bai, "GB") + size(chimeric_events, "GB") + 2) + " HDD"
#        docker: docker
#        cpu: 1
#        memory: "2GB"
#    }
#}
#

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
        set -e

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

task CreateViralFasta {
    input {
        File fasta
        String picard
        Int preemptible
        String docker
    }

    command <<<
        mv ~{fasta} virus.fasta
        # we need fasta, index, and dict in same directory
        samtools faidx virus.fasta
        java -jar ~{picard} CreateSequenceDictionary -R virus.fasta
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "1G"
        disks: "local-disk " + ceil(size(fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File viral_fasta = "virus.fasta"
        File viral_fasta_index = "virus.fasta.fai"
        File viral_dict = "virus.dict"
    }
}

task ExtractViralReads {
    input {
        File bam
        File bai
        File fasta
        String picard
        Int preemptible
        String docker
    }

    command <<<
        samtools faidx ~{fasta}

        python <<CODE

        def parse_fai(path):
            values = set()
            with open(path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        values.add(line.split('\t')[0])
            return values

        def to_txt(values, path):
            is_first = True
            with open(path, 'wt') as f:
                for val in values:
                    if not is_first:
                        f.write(' ')
                    f.write(val)
                    is_first = False

        to_txt(parse_fai('~{fasta}.fai'), 'extract.txt')
        CODE

        samtools view -b ~{bam} $(cat extract.txt) > tmp.bam
        samtools sort tmp.bam > virus.bam
        samtools index virus.bam
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "1G"
        disks: "local-disk " + ceil(size(bam, "GB")*2 + size(fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File viral_bam = "virus.bam"
        File viral_bai = "virus.bam.bai"
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
        set -e

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

#task RefineVIFOutput {
#    input {
#        File prelim_counts # InsertionSiteCandidates.abridged
#        File vif_counts # ChimericContigEvidenceAnalyzer.evidence_counts
#
#        String util_dir
#        Int preemptible
#        String docker
#    }
#
#    command <<<
#        set -e
#
#        ~{util_dir}/refine_VIF_output.Rscript \
#        --prelim_counts ~{prelim_counts} \
#        --vif_counts ~{vif_counts} \
#        --output vif.refined.tsv
#    >>>
#
#    output {
#        File refined_counts =  "vif.refined.tsv"
#    }
#
#    runtime {
#        preemptible: preemptible
#        disks: "local-disk " + ceil( size(prelim_counts, "GB")*2 + 1) + " HDD"
#        docker: docker
#        cpu: 1
#        memory: "1GB"
#    }
#}
#
#task GenomeAbundancePlot {
#    input {
#        File counts
#        String output_name
#        String title
#        String util_dir
#        Int preemptible
#        String docker
#    }
#
#    command <<<
#        set -e
#
#        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
#        --vif_report ~{counts} \
#        --title "~{title}" \
#        --output_png ~{output_name}.png
#    >>>
#
#    output {
#        File plot =  "~{output_name}.png"
#    }
#
#    runtime {
#        preemptible: preemptible
#        disks: "local-disk " + ceil( size(counts, "GB")*2 + 1) + " HDD"
#        docker: docker
#        cpu: 1
#        memory: "1GB"
#    }
#}

task VirusReport {
    input {
        File bam
        File bai
        File viral_fasta
        File insertion_site_candidates
        String util_dir
        Int preemptible
        String docker
        String memory
    }
    Int max_coverage = 100

    command <<<
        set -e

        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
        --vif_report ~{insertion_site_candidates} \
        --title "Preliminary Genome Wide Abundance" \
        --output_png vif.prelim.genome_plot.png

        # generates read_counts_summary and images
        ~{util_dir}/plot_top_virus_coverage.Rscript \
        --vif_report ~{insertion_site_candidates} \
        --bam ~{bam} \
        --output_prefix vif

        ~{util_dir}/create_insertion_site_inspector_js.py \
        --VIF_summary_tsv vif.virus_read_counts_summary.tsv \
        --json_outfile vif.virus.json

        # make bed for igvjs
        ~{util_dir}/create_virus_bed.py vif.virus_read_counts_summary.tsv vif.virus.bed

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
        File genome_abundance_plot = "vif.prelim.genome_plot.png"
        File read_counts_summary = "vif.virus_read_counts_summary.tsv"
        File read_counts_image = "vif.virus_read_counts.png"
        File read_counts_log_image = "vif.virus_read_counts_log.png"
        Array[File] virus_images = glob("vif.virus_coverage_*.png")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(bam, "GB")*2 + size(viral_fasta, "GB")*2 + 2) + " HDD"
        docker: docker
        cpu: 1
        memory: memory
    }
}

task SummaryReport {
    input {
        File prelim_counts
        File vif_counts
        File alignment_bam
        File alignment_bai
        File chim_targets_gtf
        File chim_targets_fasta
        File gtf
        Array[File] images
        String util_dir
        Int preemptible
        String docker
        String memory
    }
    Int max_coverage = 100
    String image_prefix = if(length(images)>0) then "--image " else ""
    command <<<
        set -e

        ~{util_dir}/refine_VIF_output.Rscript \
        --prelim_counts ~{prelim_counts} \
        --vif_counts ~{vif_counts} \
        --output vif.refined.tsv

        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
        --vif_report vif.refined.tsv \
        --title "Genome Wide Abundance" \
        --output_png vif.genome_plot.png

        ~{util_dir}/find_closest.py \
        -i vif.refined.tsv \
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
        --html_output vif.html

        # generate the final report
        ~{util_dir}/add_to_html.py \
        --html vif.html \
        --out vif.html \
        --image vif.genome_plot.png \
        ~{image_prefix}~{sep=' --image ' images}
    >>>

    output {
        File html = "vif.html"
        File refined_counts = "vif.refined.tsv"
        File genome_abundance_plot = "vif.genome_plot.png"
    }
    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(images, "GB") + size(alignment_bam, "GB")*2 + size(prelim_counts, "GB") + size(vif_counts, "GB") + size(chim_targets_fasta,"GB")*2 + 2) + " HDD"
        docker: docker
        cpu: 1
        memory: "16GB"
    }
}
