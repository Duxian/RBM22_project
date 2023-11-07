#######################
# 计算 pause release ratio
#######################
configfile: "config6_PRR.yaml"

peakbindTrans = config["peakbindTrans"]
bwdir = config["bwdir"]
samples = config["samples"]
outdir = config["outdir"]

transcripts = "/projects/Annotation/hg19/genecodeV19/annoFile/gencode.v19.annotation.bed12"
genes = "/projects/Annotation/hg19/genecodeV19/annoFile/Genes.bed"
gene_info = "/projects/Annotation/hg19/genecodeV19/annoFile/Genes_info.bed8"



rule all:
    input:
        expand(outdir + "Promoter_TSS-300_TSS+100.bed"),
        expand(outdir + "Genebody_TSS+300_TES-500.bed"),
        expand(outdir + "1_signal_Pro_GB/{sample}_Pro_signal.tab", sample=samples),
        expand(outdir + "2_PRR/{sample}.prr", sample=samples),
        expand(outdir + "2_PRR/" + samples[1] + ".prrFC")



# filter: gene length>2Kb; intergenic length >1Kb
rule filter_genes_len:
    input:
        trans = peakbindTrans,
    output:
        filter = outdir + "Trans.fil.bed",
        Pro = outdir + "Promoter_TSS-300_TSS+100.bed",
        GB = outdir + "Genebody_TSS+300_TES-500.bed",
    shell:
        """
        python Pyscripts/gene_filter_strand.py {input.trans} {output.filter}
        
        awk -v OFS='\\t' '{{if($6=="+") {{print $1,$2-300,$2+100,$4,$5,$6;}} else {{print $1,$3-100,$3+300,$4,$5,$6;}}}}' {output.filter} > {output.Pro}
        
        awk -v OFS='\\t' '{{if($6=="+") {{print$1,$2+300,$3-500, $4,$5,$6;}} else {{print $1,$2+500,$3-300,$4,$5,$6;}}}}' {output.filter} > {output.GB}
        """

### 计算PRR : GB/Pro
rule PRR:
    input:
        bw = bwdir + "{sample}.subtract.bw",
        Pro = outdir + "Promoter_TSS-300_TSS+100.bed",
        GB = outdir + "Genebody_TSS+300_TES-500.bed",
    output:
        Pro = outdir + "1_signal_Pro_GB/{sample}_Pro_signal.tab",
        GB = outdir + "1_signal_Pro_GB/{sample}_GB_signal.tab",
        prr = outdir + "2_PRR/{sample}.prr",
    shell:
        """
        bigWigAverageOverBed {input.bw} {input.Pro} {output.Pro}
        bigWigAverageOverBed {input.bw} {input.GB} {output.GB}
        
        python Pyscripts/prr_ratio.py {output.Pro} {output.GB} {output.prr}
        """

rule filter_PRR:
    input:
        prr = outdir + "2_PRR/{sample}.prr",
        gene = outdir + "1_signal_Pro_GB/bindPromoter_{sample}.txt",
    output:
        temp(outdir + "2_PRR/{sample}.header"),
        temp(outdir + "2_PRR/{sample}.prr.fil.tmp"),
        outdir + "2_PRR/{sample}.prr.fil",
    shell:
        """
        awk 'NR==FNR{{a[$1]=1}} NR>FNR&&a[$1]==1 {{print $0}}' {input.gene} {input.prr} > {output[1]}
        head -1 {input.prr} > {output[0]}
        cat {output[0]} {output[1]} > {output[2]}
        """


rule PRR_FC:
    input:
        KD = outdir + "2_PRR/" + samples[1] + ".prr",
        NC = outdir + "2_PRR/" + samples[0] + ".prr",
    output:
        outdir + "2_PRR/" + samples[1] + ".prrFC",
    params:
        genes = gene_info,
    shell:
        """
        python Pyscripts/prr_ratio_FC.py {input.NC} {input.KD} {params.genes} {output} 
        """

