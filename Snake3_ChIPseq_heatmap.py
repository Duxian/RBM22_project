""" get global profile heatmap for ChIP-seq """

#注意 samples 中的顺序，先KD后Control

configfile: "config3_ChIPseq_heatmap.yaml"
samples=config["SAMPLES"]
bwdir = config["bwdir"]
outdir = config["outdir"]
name = config["name"]
genes = config["genes"]
groups = config["groups"]

fa = "/projects/Genome/hg19/hg19.fasta"

rule all:
    input:
        outdir + "/" + name + "/0_substract_KD-NC/substract_KD-NC.bw",
        expand(outdir + "/" + name + "/{group}/1_makeTab/" + name + "_KD-NC-sub.gz", group=groups),
        expand(outdir + "/" + name + "/{group}/2_plotHeatmap/" + name + "_KD-NC-sub.pdf", group=groups),
        expand(outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.gz", group=groups),
        expand(outdir + "/{group}/2_plotHeatmap_aroundTES/" + name + "_KD-NC-sub.pdf", group=groups),
        expand(outdir + "/{group}/3_countSignal_aroundTES/" + name + "_sub.average.txt",group=groups),

## KD - NC
rule substract_kd_ctl:
    input:
        kdbw = bwdir + samples[1] + ".subtract.bw",
        ctlbw = bwdir + samples[0] + ".subtract.bw",
    output:
        subbw = outdir + "/" + name + "/0_substract_KD-NC/substract_KD-NC.bw",
    threads: 40
    shell:
        """
        bigwigCompare -b1 {input.kdbw} -b2 {input.ctlbw}  -bs 1 --skipZeroOverZero --operation subtract -p {threads} -o {output.subbw} -of "bigwig"
        """


############################################################## for heatmap on genes

rule computeMatrix:
    input:
        bw = expand(bwdir + "{sample}.subtract.bw",sample=samples),
        sub1 = outdir + "/" + name + "/0_substract_KD-NC/substract_KD-NC.bw",
        gene = lambda wildcards: genes[groups.index(wildcards.group)],
    output:
        tab = outdir + "/" + name + "/{group}/1_makeTab/" + name + "_KD-NC-sub.tab",
        gz = outdir + "/" + name + "/{group}/1_makeTab/" + name + "_KD-NC-sub.gz",
        sortReg = outdir + "/" + name + "/{group}/1_makeTab/" + name + "_KD-NC-sub.sortRegion.bed",
    threads: 40
    params:
        b = 2000,
        a = 0,
        m = 4000,
    shell:
        """
        computeMatrix scale-regions -S {input.bw} {input.sub1} -R {input.gene} \
        --samplesLabel "NC" "KD" "KD-NC" -b {params.b} -a {params.a} -m {params.m} \
        --sortRegions descend --missingDataAsZero --skipZeros -p {threads} \
        -o {output.gz} --outFileNameMatrix {output.tab} --outFileSortedRegions {output.sortReg} 
        """



rule plotheatmap:
    input:
        gzAll = outdir  + "/" + name + "/{group}/1_makeTab/" + name + "_KD-NC-sub.gz",
    output:
        pdf = outdir  + "/" + name + "/{group}/2_plotHeatmap/" + name + "_KD-NC-sub.pdf",
        tab = outdir  + "/" + name + "/{group}/2_plotHeatmap/" + name + "_heatmap.matrix.gz",
    shell:
        """
        plotHeatmap -m {input.gzAll} -o {output.pdf} \
        --sortRegions keep \
        --colorMap  Reds Blues Greys RdBu RdGy \
        --whatToShow 'plot, heatmap and colorbar' \
        --missingDataColor 'white' \
        --dpi 300 \
        --zMax 3 3 3 0.8 0.8 --zMin 0 0 0 -0.8 -0.8\
        --outFileNameMatrix {output.tab} \
        --yAxisLabel "reads density" \
        --plotFileFormat "pdf" \
        --heatmapHeight 10
        """




##########################################################  for heatmap on TES
# [TES-2Kb, TES+10Kb]
rule computeMatrix_TES:
    input:
        bw = expand(bwdir + "{sample}.bw",sample=samples),
        sub = outdir + "/0_substract_KD-NC/substract_KD-NC.bw",
        gene = lambda wildcards: genes[groups.index(wildcards.group)],
    output:
        tab = outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.tab",
        gz = outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.gz",
        sortReg = outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.sortRegion.bed",
    threads: 60
    params:
        b = 2000,
        a = 10000,
    shell:
        """
        computeMatrix reference-point \
        -S {input.bw} {input.sub} \
        -R {input.gene} \
        --samplesLabel "KD" "NC" "KD-NC" \
        -b {params.b} -a {params.a}  --referencePoint TES \
        --sortRegions descend --missingDataAsZero --skipZeros \
        -p {threads} \
        -o {output.gz} \
        --outFileNameMatrix {output.tab} --outFileSortedRegions {output.sortReg} 
        """

rule plotheatmap_TES:
    input:
        gzAll = outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.gz",
    output:
        pdf = outdir + "/{group}/2_plotHeatmap_aroundTES/" + name + "_KD-NC-sub.pdf",
        tab = outdir + "/{group}/2_plotHeatmap_aroundTES/" + name + "_heatmap.matrix.gz",
    shell:
        """
        plotHeatmap -m {input.gzAll} -o {output.pdf} \
        --sortRegions keep \
        --colorMap Oranges Purples PuOr_r \
        --whatToShow 'plot, heatmap and colorbar' \
        --missingDataColor 'white' \
        --dpi 300 \
        --zMax 0.04 0.06 0.04 \
        --zMin 0 0 -0.04 \
        --outFileNameMatrix {output.tab} \
        --yAxisLabel "reads density" \
        --plotFileFormat "pdf" \
        --heatmapHeight 10
        """



# make signal for each position
rule count_on_position:
    input:
        tab = outdir + "/{group}/1_makeTab_aroundTES/" + name + "_KD-NC-sub.tab",
    output:
        scale = temp(outdir + "/{group}/3_countSignal_aroundTES/" + name + "_KD-NC-sub.tab"),
        kd = outdir + "/{group}/3_countSignal_aroundTES/" + name + "_KD.average.txt",
        ctl = outdir + "/{group}/3_countSignal_aroundTES/" + name + "_NC.average.txt",
        sub = outdir + "/{group}/3_countSignal_aroundTES/" + name + "_sub.average.txt",
    params: 3
    shell:
        """
        sed '1,3d' {input.tab} > {output.scale}
        python scripts/countAverage_threeSample.py {output.scale} {params} {output.kd} {output.ctl} {output.sub}
        """


