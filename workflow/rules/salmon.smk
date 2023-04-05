localrules: copy_salmon_indices_to_mount

rule copy_salmon_indices_to_mount:
    """
    Copy to a dir that singularity will see.
    """
    input:
        index = refs_wf("results/salmon-idx/index/"),
        aux = refs_wf("results/salmon-idx/transcripts_and_consensus_tes.aux.txt"),
        fasta = refs_wf("results/combined-anno/transcripts-plus-tes.fasta.gz"),
        gtf = refs_wf("results/combined-anno/transcripts-plus-tes.gtf"),
        tx2gene = refs_wf("results/combined-anno/transcripts-plus-tes.tx2txsymbol.tsv"),
    output:
        index = directory("results/indices/vanilla_salmon_tes_transcripts/index/"),
        aux = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.aux.txt",
        fasta = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.fasta.gz",
        gtf = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf",
        tx2gene = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2txsymbol.tsv",
    priority: 51
    shell:
        """
        cp -r {input.index} {output.index}
        cp {input.aux} {output.aux}
        cp {input.fasta} {output.fasta}
        cp {input.gtf} {output.gtf}
        cp {input.tx2gene} {output.tx2gene}
        """

rule salmon_quant_se_vanilla:
    """
    We ran salmon v1.5.2 quant directly on trimmed fastqs with options
    "<replace with options actually used>".
    """
    input:
        reads = lambda wc: [rules.fastp_trim_pe.output.r1,rules.fastp_trim_pe.output.r2] if is_paired_end(wc.sample) else rules.fastp_trim_se.output.r1,
        idx = rules.copy_salmon_indices_to_mount.output.index,
        aux = rules.copy_salmon_indices_to_mount.output.aux,
    output:
        sf = "results/quantification/vanilla_salmon_tes_transcripts/quant/{sample}/{quant_rep}/quant.sf",
        #sam = temp("results/quantification/vanilla_salmon_tes_transcripts/{sample}/{quant_rep}/alignments/alignments.sam")
    resources:
        time=60,
        mem=20000,
        cpus=8
    params:
        fq_args = lambda wc: " ".join([" ".join(x) for x in zip(["-1","-2"],expand("results/reads/fastp-pe/{s}_r{e}.trimmed.fq.gz",s=wc.sample,e=[1,2]))]) if is_paired_end(wc.sample) else "-r results/reads/fastp-se/{s}_r1.trimmed.fq.gz".format(s=wc.sample),
        libtype = config.get("SALMON_VANILLA_LIBTYPE"),
        bootstraps = "--numBootstraps {n}".format(n=config.get("SALMON_VANILLA_BOOTSTRAPS")) if config.get("SALMON_VANILLA_BOOTSTRAPS") else "",
        seqbias = "--seqBias" if config.get("SALMON_VANILLA_SEQBIAS") else "",
        gcbias = "--gcBias" if config.get("SALMON_VANILLA_GCBIAS") else "",
        posbias = "--posBias" if config.get("SALMON_VANILLA_POSBIAS") else "",
        auxtargetfile = "--auxTargetFile {x}".format(x="results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.aux.txt") if config.get("SALMON_VANILLA_USE_AUXTARGETS") else "",
        softclip = "--softclip" if config.get("SALMON_VANILLA_SOFTCLIP") else "",
        incompat = "--incompatPrior 0.0" if config.get("SALMON_VANILLA_INCOMPATPRIOR_0") else ""
    threads:
        8
    log:
        "results/logs/salmon_quant_se_vanilla/{sample}.{quant_rep}.txt"
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    priority: 51
    shell:
        """
        salmon quant --index {input.idx} \
            --libType {params.libtype} \
            {params.fq_args} \
            {params.bootstraps} \
            {params.seqbias} \
            {params.gcbias} \
            {params.posbias} \
            {params.auxtargetfile} \
            {params.softclip} \
            {params.incompat} \
            --writeUnmappedNames \
            -d \
            -p {threads} \
            -o $(dirname {output.sf})/ 2> {log}
        """

rule vanilla_salmon_linked_txome:
    input:
        idx = rules.copy_salmon_indices_to_mount.output.index,
        fasta = rules.copy_salmon_indices_to_mount.output.fasta,
        gtf = rules.copy_salmon_indices_to_mount.output.gtf,
        tx2gene = rules.copy_salmon_indices_to_mount.output.tx2gene,
    output:
        json = "results/quantification/vanilla_salmon_tes_transcripts/tximeta.json",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-tximeta:1.10.0--r41hdfd78af_0"
    resources:
        time=10,
        mem=4000,
        cpus=1
    params:
        flybase_release = config.get("FLYBASE_RELEASE"),
        genome_version = config.get("GENOME_VERSION"),
    log:
        "results/logs/vanilla_salmon_linked_txome/log.txt"
    script:
        "../scripts/vanilla_salmon_linked_txome.R"

rule vanilla_salmon_tximeta:
    """
    tximeta was used to summarize expression estimates and include metadata in summarizedExperiment objects.
    """
    input:
        json = rules.vanilla_salmon_linked_txome.output.json,
        samples = "config/sample_table.csv", # TODO make this the output of the anno rule
        salmon_files = expand("results/quantification/vanilla_salmon_tes_transcripts/quant/{s}/{{quant_rep}}/quant.sf", s=SAMPLES),
        tx2gene = rules.copy_salmon_indices_to_mount.output.tx2gene,
        pipeline_meta = rules.get_pipeline_info.output,
        # not manually used in this script, but must still exist
        idx = rules.copy_salmon_indices_to_mount.output.index,
        fasta = rules.copy_salmon_indices_to_mount.output.fasta,
        gtf = rules.copy_salmon_indices_to_mount.output.gtf,
    output:
        salmon = "results/quantification/vanilla_salmon_tes_transcripts/se.gene.{quant_rep}.rds",
        salmon_tx ="results/quantification/vanilla_salmon_tes_transcripts/se.transcript.{quant_rep}.rds",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-tximeta:1.10.0--r41hdfd78af_0"
    params:
        counts_from_abundance_salmon = config.get("SALMON_VANILLA_TXIMPORT_COUNTSFROMABUNDANCE"),
        counts_from_abundance_salmon_txout = config.get("SALMON_VANILLA_TXIMPORT_COUNTSFROMABUNDANCE_TXOUT"),
    resources:
        time=60,
        mem=128000,
        cpus=1
    log:
        "results/logs/vanilla_salmon_tximeta/log.{quant_rep}.txt"
    script:
        "../scripts/vanilla_salmon_tximeta.R"