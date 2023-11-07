"""make track from bam file for ChIP-seq"""

# source activate py3

configfile: "config2_makeTrack.yaml"     
SAMPLES = config["SAMPLES"]
Inputs = config["Inputs"]
Bamdir = config["Bamdir"]
outdir = config["outdir"]

hg19_genomeSize = "projects/Genome/hg19/hg19.chrom.sizes"


rule all:
    input:
        expand(outdir + "/1_bw_nor_CPM/{sample}.nor.bw", sample=SAMPLES + Inputs),
        expand(outdir + "/2_bw_subtract_IP-Input/{sample}.subtract.bw", sample=SAMPLES),
        expand(outdir + "/2_bw_subtract_IP-Input_adjust/{sample}.subtract.bw", sample=SAMPLES),
        expand(outdir + "/3_bed/{sample}.bed", sample=SAMPLES + Inputs),




#1. normalize bigwig file from bam  
rule bam2bw:
    input:
        bam = Bamdir + "{sample}.nodup.bam",    
    output:
        bw = outdir + "/1_bw_nor_CPM/{sample}.nor.bw",
    threads: 60,
    shell:
        """
        bamCoverage --bam {input.bam} -o {output.bw} -p {threads} --binSize 1 --normalizeUsing CPM 
        """

#2. subtract for IP -Input
rule subtract_IP_Input:
    input:
        IP = outdir + "/1_bw_nor_CPM/{sample}.nor.bw",
    params:
        IgG = lambda wildcards: [outdir + "/1_bw_nor_CPM/" + Inputs[SAMPLES.index(wildcards.sample)] + ".nor.bw"],
    threads: 60,
    output:
        bw = outdir + "/2_bw_subtract_IP-Input/{sample}.subtract.bw",
    shell:
        """
        bigwigCompare -b1 {input.IP} -b2 {params.IgG} --skipNonCoveredRegions --operation subtract --binSize 1 -p {threads} \
        -o {output.bw} --skipNonCoveredRegions
        """

#convert negative value to 0
rule convert_value:
    input:
        bw = outdir + "/2_bw_subtract_IP-Input/{sample}.subtract.bw",
    output:
        bg = outdir + "/2_bw_subtract_IP-Input_adjust/{sample}.subtract.bg",
        newbg = temp(outdir + "/2_bw_subtract_IP-Input_adjust/{sample}.subtract.newbg"),
        sortbg = temp(outdir + "/2_bw_subtract_IP-Input_adjust/{sample}.subtract.sort.bg"),
        bw = outdir + "/2_bw_subtract_IP-Input_adjust/{sample}.subtract.bw",
    params:
        genomeSize = hg19_genomeSize,
    shell:
        """
        bigWigToBedGraph {input.bw} {output.bg}   && \
        awk '{{if($4<0) {{print $1"\\t"$2"\\t"$3"\\t""0";}} else {{print $1"\\t"$2"\\t"$3"\\t"$4;}}}}' {output.bg} > {output.newbg}   && \
        LC_COLLATE=C sort -k1,1 -k2n {output.newbg} > {output.sortbg}  && \
        bedGraphToBigWig {output.sortbg} {params.genomeSize} {output.bw}
        """

#3.bam to bed
rule bam2bed:
    input:
        bam = Bamdir + "{sample}.sbam",
    output:
        bed = outdir + "/3_bed/{sample}.bed",
    shell:
        """
        bedtools bamtobed -split -bed12 -i {input.bam} > {output.bed}
        """


