""" compute readthrough index"""

configfile: "config7_RT.yaml"

samples= config["samples"]
bwdir = config["bwdir"]
DS = config["DS"]
outdir = config["outdir"]

hg19_genomeSize = "/projects/Elongation_RBM22_project/data/hg19Annotation/hg19.clean.chrom.sizes"
hg19_Genes_info = "/projects/Elongation_RBM22_project/data/hg19Annotation/Genes_info.bed8"
hg19_gff3 = "/projects/Annotation/hg19/genecodeV19/gencode.v19.annotation.gff3"
PC="/projects/Annotation/hg19/genecodeV19/annoFile/PC_genes_info.bed8"



rule all:
    input:
        outdir + "lastExon_eachGenes.txt",
        outdir + "PC_lastExon.rev.txt",
        expand(outdir + "1_signal_DS_LastE/{sample}_DS_signal.fwd.tab", sample=samples),
        expand(outdir + "2_RI/{sample}.RI", sample=samples),



# 寻找每个基因最远端的last exon
rule get_lastExon:
    input:
        hg19_gff3,
    output:
        outdir + "lastExon_eachGenes.txt",
    shell:
        """
        python Pyscripts/getLastExon.py  {input} {output} 
        """

# 筛选protein-coding gene的last exon 和 downstream intergenic region对应
rule PC_match:
    input:
        DS = DS,
        PC = PC,
        LastE = outdir + "lastExon_eachGenes.txt",
    output:
        temp(outdir +  "gene_downstream_region.PC.txt.tmp"),
        outdir +  "PC_downstream_region.fwd.txt",
        outdir + "PC_downstream_region.rev.txt",
        outdir + "PC_lastExon.fwd.txt",
        outdir + "PC_lastExon.rev.txt",
    shell:
        """
        sed 's/_inter//g' {input.DS} |awk -v OFS='\\t' '{{print $2,$3,$4,$1,"0",$5}}' > {output[0]}
        
        awk 'NR==FNR{{a[$4]=1}} NR>FNR&&a[$4]==1 {{print $0}}' {input.PC} {output[0]} |awk '$6=="+"{{print}}' > {output[1]}
        awk 'NR==FNR{{a[$4]=1}} NR>FNR&&a[$4]==1 {{print $0}}' {input.PC} {output[0]} |awk '$6=="-"{{print}}' > {output[2]}
        
        awk 'NR==FNR{{a[$4]=1}} NR>FNR&&a[$7]==1 {{print $0}}' {input.PC} {input.LastE} | awk -v OFS='\\t' '$5=="+"{{print $1,$2,$3,$7,"0",$5}}' > {output[3]}
        awk 'NR==FNR{{a[$4]=1}} NR>FNR&&a[$7]==1 {{print $0}}' {input.PC} {input.LastE} | awk -v OFS='\\t' '$5=="-"{{print $1,$2,$3,$7,"0",$5}}' > {output[4]}
        """


rule count_density:
    input:
        fwd = bwdir + "{sample}_forward.bw",
        rev = bwdir + "{sample}_reverse.bw",
	pc_DS_f = outdir + "PC_downstream_region.fwd.txt",
        pc_DS_r = outdir + "PC_downstream_region.rev.txt",
        pc_E_f = outdir + "PC_lastExon.fwd.txt",
        pc_E_r = outdir + "PC_lastExon.rev.txt",
    output:
        DS_fwd = outdir + "1_signal_DS_LastE/{sample}_DS_signal.fwd.tab",
        DS_rev = outdir + "1_signal_DS_LastE/{sample}_DS_signal.rev.tab",
        LastE_fwd = outdir + "1_signal_DS_LastE/{sample}_LastExon_signal.fwd.tab",
        LastE_rev = outdir + "1_signal_DS_LastE/{sample}_LastExon_signal.rev.tab",
    shell:
        """
        bigWigAverageOverBed {input.fwd} {input.pc_DS_f} {output.DS_fwd}
        bigWigAverageOverBed {input.rev} {input.pc_DS_r} {output.DS_rev}
        bigWigAverageOverBed {input.fwd} {input.pc_E_f} {output.LastE_fwd}
        bigWigAverageOverBed {input.rev} {input.pc_E_r} {output.LastE_rev}
        """

rule RI:
    input:
        DS_fwd = outdir + "1_signal_DS_LastE/{sample}_DS_signal.fwd.tab",
        DS_rev = outdir + "1_signal_DS_LastE/{sample}_DS_signal.rev.tab",
        LastE_fwd = outdir + "1_signal_DS_LastE/{sample}_LastExon_signal.fwd.tab",
        LastE_rev = outdir + "1_signal_DS_LastE/{sample}_LastExon_signal.rev.tab",
    output:
        outdir + "2_RI/{sample}_DS_signal.tab",
        outdir + "2_RI/{sample}_LastExon_signal.tab",
        outdir + "2_RI/{sample}.RI",
    shell:
        """
        cat {input.DS_fwd} {input.DS_rev} > {output[0]}
        cat {input.LastE_fwd} {input.LastE_rev} > {output[1]}
        
        python Pyscripts/RI_ratio.py {output[0]} {output[1]} {output[2]} 
        """




