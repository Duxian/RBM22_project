"""call peak by MACS2 and find precise transcript TSS bound by RNAPII"""

configfile: "config4.yaml"

chip = config["chip"]
Input = config["Input"]
names = config["names"]
Bamdir = config["Bamdir"]
outdir = config["outdir"]



genome = "/projects/Genome/hg19/hg19.chrom.sizes"
genes = "/projects/Annotation/hg19/genecodeV19/annoFile/Genes.bed"
transcripts = "/projects/Annotation/hg19/genecodeV19/annoFile/gencode.v19.annotation.bed12"
ID = "/projects/Annotation/hg19/genecodeV19/annoFile/trans2gene_ID.txt"


rule all:
    input:
        expand(outdir + "/1_MACS2/narrowPeak/{name}_peaks.narrowPeak", name=names),
        expand(outdir + "/1_MACS2/broadPeak/{name}_peaks.broadPeak", name=names),
        expand(outdir + "/1_MACS2/narrowPeak/{name}_peaks.filt.narrowPeak.bed", name=names),
        expand(outdir + "/1_MACS2/narrowPeak/{name}_summits.fil.bed", name=names),
        expand(outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet", name=names),
        expand(outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet.addID", name=names),
        expand(outdir + "/2_closetTSS/narrowPeak/{name}_preciseTrans.bed12", name=names)



####PART1: call peaks

#1.call peak for POLR2G chip of KD RBM22 and NC
rule MACS2_POLR2GChip:
    input:
        chip = lambda wildcards: (Bamdir + chip[names.index(wildcards.name)] + ".nodup.bam"),
        Input = lambda wildcards: (Bamdir + Input[names.index(wildcards.name)] + ".nodup.bam"),
    output:
        outdir + "/1_MACS2/narrowPeak/{name}_peaks.narrowPeak",
        outdir + "/1_MACS2/broadPeak/{name}_peaks.broadPeak",
        outdir + "/1_MACS2/narrowPeak/{name}_summits.bed",
    params:
        prefix = "{name}",
        narrow_out = directory(outdir + "/1_MACS2/narrowPeak/"),
        broad_out = directory(outdir + "/1_MACS2/broadPeak/"),
    log:
        log1 = outdir + "/1_MACS2/narrowPeak/{name}.narrlog",
        log2 = outdir + "/1_MACS2/broadPeak/{name}.broadlog",
    shell:
        """
        macs2 callpeak -t {input.chip} -c {input.Input} \
        -f BAM -g hs --nomodel --keep-dup all -n {params.prefix} \
        --max-gap 50 --outdir {params.narrow_out}  1>{log.log1} 2>&1  && \
        
        macs2 callpeak -t {input.chip} -c {input.Input} \
        -f BAM -g hs --nomodel --keep-dup all -n {params.prefix} \
        --max-gap 50 --outdir {params.broad_out}  --broad 1>{log.log2} 2>&1
        """

# filter for more significant peak : FC>4 and q-value<0.00001
rule filter_peak_POLR2GChip:
    input:
        outdir + "/1_MACS2/narrowPeak/{name}_peaks.narrowPeak",
    output:
        outdir + "/1_MACS2/narrowPeak/{name}_peaks.filt.narrowPeak",
        outdir + "/1_MACS2/narrowPeak/{name}_peaks.filt.narrowPeak.bed",
    shell:
        """
        awk '$7>=4&&$9>=5{{print}}' {input} > {output[0]}
        cut -f 1-6 {output[0]} > {output[1]}
        """

rule filter_for_summits:
    input:
        outdir + "/1_MACS2/narrowPeak/{name}_peaks.filt.narrowPeak",
        outdir + "/1_MACS2/narrowPeak/{name}_summits.bed",
    output:
        outdir + "/1_MACS2/narrowPeak/{name}_summits.fil.bed",
    shell:
        """
        awk 'NR==FNR{{a[$4]=1}} NR>FNR&&a[$4]==1 {{print $0}}' {input[0]} {input[1]} > {output}
        """



###PART2 get closet TSS according peak : 区间范围越窄，越能找到准确的TSS

rule closet:
    input:
        summit = outdir + "/1_MACS2/narrowPeak/{name}_summits.fil.bed",
    output:
        temp(outdir + "/2_closetTSS/narrowPeak/{name}_summits.extd400.bed"),
        outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet",
    params:
        hg19 = genome,
        transcripts = transcripts,
        TSS = outdir + "/2_closetTSS/narrowPeak/Transcripts_TSS_add100.bed",
    shell:
        """
        awk '{{if($6=="+") {{print $1"\\t"$2"\\t"$2+100"\\t"$4"\\t"$5"\\t"$6;}} \
        else {{print $1"\\t"$3-100"\\t"$3"\\t"$4"\\t"$5"\\t"$6;}}}}'  {params.transcripts} |sort -k1,1 -k2n -k3n -k6,6 -u > {params.TSS}
        
        bedtools slop -b 100 -i {input.summit} -g {params.hg19} |sort -k1,1 -k2n > {output[0]}
        bedtools closest -a {output[0]} -b {params.TSS} -d > {output[1]}
        """


rule add_GeneID:
    input:
        outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet",
    output:
        temp(outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet.tmp"),
        outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet.addID",
    params:
        ID = ID,
        sortID = outdir + "/2_closetTSS/narrowPeak/trans2gene_ID.sort.txt",
    shell:
        """
        sort -k9,9 {input} > {output[0]}
        sort -k1,1 {params.ID} > {params.sortID}
        join -1 9 -2 1 -t $'\\t' {output[0]} {params.sortID} > {output[1]}
        """


# find precise transcript TSS for each gene that P2_NC bind
rule find_preciseTSS:
    input:
        outdir + "/2_closetTSS/narrowPeak/{name}_Peaks_TranscriptsTSS.closet.addID",
    output:
        fwd = temp(outdir + "/2_closetTSS/narrowPeak/{name}_preciseTSS.Fwd.txt"),
        rev = temp(outdir + "/2_closetTSS/narrowPeak/{name}_preciseTSS.Rev.txt"),
        cat = outdir + "/2_closetTSS/narrowPeak/{name}_preciseTSS.ID",
        bed = outdir + "/2_closetTSS/narrowPeak/{name}_preciseTrans.bed12",
    params:
        transcripts = transcripts,
    shell:
        """
        awk '$12<200{{print}}' {input} | awk '$11=="+"{{print}}' |sort -k13,13 -k12n -k8n |sort -k13 -u > {output.fwd}
        awk '$12<200{{print}}' {input} | awk '$11=="-"{{print}}' |sort -k13,13 -k12n -k9nr |sort -k13 -u > {output.rev}
        
        cat {output.fwd} {output.rev} |cut -f 1,13,14,15 > {output.cat}
        awk 'NR==FNR{{a[$1]=1}} NR>FNR&&a[$4]==1 {{print $0}}' {output.cat} {params.transcripts} > {output.bed}
        """


