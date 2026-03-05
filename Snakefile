configfile: "config.yaml"

IN_DIR = config["input_dir"]
AN_DIR = config["analysis_dir"]

# Auto-detect samples
SAMPLES = glob_wildcards(f"{IN_DIR}" + "/{sample}.fastq.gz").sample

rule all:
    conda: 
        "envs/denovo.yaml"
    input:
        expand(f"{AN_DIR}/variants_over1pct/{{sample}}_variants.consolidated.fasta", sample=SAMPLES),
        expand(f"{AN_DIR}/variants_over1pct/{{sample}}_variants.consolidated.abundance.tsv", sample=SAMPLES),
        expand(f"{AN_DIR}/variants_over1pct/{{sample}}_read_nseqs_summary.tsv", sample=SAMPLES)

rule filter_chopper:
    conda:
        "envs/denovo.yaml"
    input:
        fastq = f"{IN_DIR}/{{sample}}.fastq.gz"
    output:
        filt = temp(f"{AN_DIR}/filt/{{sample}}_filtered.fastq.gz")
    params:
        q = config["quality"],
        minlen = config["minlen"],
        maxlen = config["maxlen"]
    shell:
        r"""
        mkdir -p {AN_DIR}/filt
        chopper -q {params.q} -l {params.minlen} --maxlength {params.maxlen} -i {input.fastq} \
        | gzip > {output.filt}
        """

rule nseq_filtered_reads:
    input:
        filt = f"{AN_DIR}/filt/{{sample}}_filtered.fastq.gz"
    output:
        total = temp(f"{AN_DIR}/filt/{{sample}}_filtered.total.txt")
    shell:
        r"""
        zcat < {input.filt} | awk 'NR%4==1{{n++}}END{{print n}}' > {output.total}
        """

rule cluster_vsearch:
    conda:
        "envs/denovo.yaml"
    input:
        filt = f"{AN_DIR}/filt/{{sample}}_filtered.fastq.gz"
    output:
        centroids = temp(f"{AN_DIR}/centroids/{{sample}}_centroids.99.fasta"),
        consensus = temp(f"{AN_DIR}/consensus/{{sample}}_consensus.99.fasta")
    params:
        id = config["percid"],
        threads = config["threads"]["vsearch"]
    shell:
        r"""
        mkdir -p {AN_DIR}/centroids {AN_DIR}/consensus
        vsearch --cluster_fast {input.filt} \
          --id {params.id} \
          --centroids {output.centroids} \
          --consout {output.consensus} \
          --sizeout \
          --threads {params.threads}
        """

rule filter_consensus_by_abundance:
    input:
        consensus = f"{AN_DIR}/consensus/{{sample}}_consensus.99.fasta",
        total = f"{AN_DIR}/filt/{{sample}}_filtered.total.txt"
    output:
        over = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta")
    params:
        threshold_frac = config["threshold_frac"]
    shell:
        r"""
        mkdir -p {AN_DIR}/variants_over1pct
        TOTAL=$(cat {input.total})
        THRESH=$(awk -v t="$TOTAL" -v thr="{params.threshold_frac}" 'BEGIN{{printf "%d\n", (thr*t)}}')

        awk -v thr="$THRESH" '
        /^>/ {{
            keep=0
            split($0,a,"size=")
            if(length(a)>1) {{
                split(a[2],b,";")
                size=b[1]
                if(size+0 >= thr) keep=1
            }}
        }}
        {{ if(keep) print }}
        ' {input.consensus} > {output.over}
        """

rule abundance_table_over_threshold:
    input:
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta",
        total = f"{AN_DIR}/filt/{{sample}}_filtered.total.txt"
    output:
        tsv = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.abundance.tsv")
    shell:
        r"""
        TOTAL=$(cat {input.total})
        awk -v total="$TOTAL" '
        /^>/ {{
          name=$0; sub(/^>/,"",name)
          size=0
          split($0,a,"size=")
          if(length(a)>1) {{
              split(a[2],b,";")
              size=b[1]+0
          }}
          if(size>0) printf "%s\t%d\t%.6f\n", name, size, size/total
        }}' {input.over} | sort -k2,2nr > {output.tsv}
        """

checkpoint nseq_over_threshold:
    input:
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta"
    output:
        nseq = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_nseq.txt")
    shell:
        r"""
        grep -c '^>' {input.over} > {output.nseq} || echo 0 > {output.nseq}
        """

def needs_full_path(wc):
    ck = checkpoints.nseq_over_threshold.get(sample=wc.sample)
    with open(ck.output.nseq) as f:
        c = int(f.read().strip() or "0")
    return c >= 2

rule align_mafft:
    input:
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta"
    output:
        aln = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.aln.fasta")
    log:
        f"{AN_DIR}/logs/{{sample}}.mafft.log"
    conda:
        "envs/denovo.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {AN_DIR}/logs

        N=$(grep -c '^>' {input.over} || true)

        if [ "$N" -lt 2 ]; then
            # Nothing to align; pass-through so downstream can still work
            cp {input.over} {output.aln}
            echo "SKIP: <2 sequences (N=$N). Input copied to alignment output." > {log}
        else
            mafft --preservecase --auto {input.over} > {output.aln} 2> {log}
        fi
        """
rule endtrim_alignment:
    input:
        aln = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.aln.fasta",
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta"  # to count sequences robustly
    output:
        endtrim = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.endtrim.aln.fasta")
    log:
        f"{AN_DIR}/logs/{{sample}}.endtrim.log"
    conda:
        "envs/denovo.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {AN_DIR}/logs

        N=$(grep -c '^>' {input.over} || true)

        if [ "$N" -lt 2 ]; then
            # Nothing to trim; pass through the "alignment" as-is
            cp {input.aln} {output.endtrim}
            echo "SKIP: <2 sequences (N=$N). Alignment copied to endtrim output." > {log}
        else
            python3 scripts/endtrim_alignment.py {input.aln} {output.endtrim} > {log} 2>&1
        fi
        """

rule ungap_endtrimmed:
    input:
        endtrim = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.endtrim.aln.fasta",
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta"
    output:
        ungapped = temp(f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.endtrim.ungapped.fasta")
    log:
        f"{AN_DIR}/logs/{{sample}}.ungap.log"
    conda:
        "envs/denovo.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {AN_DIR}/logs

        # how many sequences exist in the true pre-alignment file?
        N=$(grep -c '^>' {input.over} || true)

        # if there are <2 sequences, there is no meaningful alignment/trim/ungap workflow
        if [ "$N" -lt 2 ]; then
            # Create a valid fasta output:
            # - if endtrim has content, copy it
            # - else, copy the original over-threshold fasta (already valid)
            if [ -s {input.endtrim} ]; then
                cp {input.endtrim} {output.ungapped}
                echo "SKIP: <2 sequences (N=$N). Copied endtrim to ungapped." > {log}
            else
                cp {input.over} {output.ungapped}
                echo "SKIP: <2 sequences (N=$N). Endtrim empty; copied original over-threshold fasta." > {log}
            fi
        else
            python3 scripts/ungap_fasta.py {input.endtrim} {output.ungapped} > {log} 2>&1
        fi
        """

rule make_consolidated_variants:
    input:
        nseq = f"{AN_DIR}/variants_over1pct/{{sample}}_nseq.txt",
        over = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.fasta",
        over_tsv = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.abundance.tsv",
        ungapped = f"{AN_DIR}/variants_over1pct/{{sample}}_consensus_over1pct.endtrim.ungapped.fasta",
        total = f"{AN_DIR}/filt/{{sample}}_filtered.total.txt"
    output:
        fasta = f"{AN_DIR}/variants_over1pct/{{sample}}_variants.consolidated.fasta",
        tsv = f"{AN_DIR}/variants_over1pct/{{sample}}_variants.consolidated.abundance.tsv"
    log:
        f"{AN_DIR}/logs/{{sample}}.consolidate.log"
    conda:
        "envs/denovo.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {AN_DIR}/logs

        # Read nseq robustly (strip whitespace, default 0)
        nseq=$(tr -d '[:space:]' < {input.nseq} || true)
        if [ -z "$nseq" ]; then nseq=0; fi

        echo "nseq=$nseq" > {log}
        echo "over_bytes=$(wc -c < {input.over})" >> {log}
        echo "ungapped_bytes=$(wc -c < {input.ungapped} 2>/dev/null || echo 0)" >> {log}
        echo "vsearch=$(command -v vsearch || echo NOT_FOUND)" >> {log}

        if [ "$nseq" -lt 2 ]; then
            cp {input.over} {output.fasta}
            cp {input.over_tsv} {output.tsv}
            echo "SKIP: nseq<2, copied over-threshold outputs." >> {log}
            exit 0
        fi

        # If ungapped is empty for some reason, fall back to over-threshold fasta
        if [ ! -s {input.ungapped} ]; then
            echo "WARN: ungapped fasta empty; falling back to over-threshold fasta for derep." >> {log}
            src={input.over}
        else
            src={input.ungapped}
        fi

        # Dereplicate and sum sizes
        vsearch --derep_fulllength "$src" \
          --output {output.fasta} \
          --sizein --sizeout >> {log} 2>&1

        TOTAL=$(cat {input.total})
        awk -v total="$TOTAL" '
        /^>/ {{
          name=$0; sub(/^>/,"",name)
          size=0
          split($0,a,"size=")
          if(length(a)>1) {{ split(a[2],b,";"); size=b[1]+0 }}
          if(size>0) printf "%s\t%d\t%.6f\n", name, size, size/total
        }}' {output.fasta} | sort -k2,2nr > {output.tsv}

        echo "OK: consolidated + abundance table written." >> {log}
        """

rule read_nseqs_summary:
    input:
        consolidated = f"{AN_DIR}/variants_over1pct/{{sample}}_variants.consolidated.fasta",
        total = f"{AN_DIR}/filt/{{sample}}_filtered.total.txt"
    output:
        tsv = f"{AN_DIR}/variants_over1pct/{{sample}}_read_nseqs_summary.tsv"
    shell:
        r"""
        TOTAL=$(cat {input.total})
        CLUSTERED=$(awk '
        /^>/ {{
          split($0,a,"size=")
          if(length(a)>1) {{ split(a[2],b,";"); s+=b[1]+0 }}
        }}
        END{{ if(s=="") s=0; print s }}
        ' {input.consolidated})
        NON_CLUSTERED=$((TOTAL - CLUSTERED))
        {{ echo -e "TOTAL\tCLUSTERED\tNON_CLUSTERED"; echo -e "${{TOTAL}}\t${{CLUSTERED}}\t${{NON_CLUSTERED}}"; }} > {output.tsv}
        """