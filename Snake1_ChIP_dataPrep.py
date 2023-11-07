"""
# data process for  ChIP-seq data
"""

# source activate python3 ; then source deactivate

configfile: "config1_ChIPseq_dataPrep.yaml"
samples = config["SAMPLES"]
DataDir = config["datadir"]
outdir = config["outdir"]

idx = "/projects/Genome/bt2_idx/hg19/hg19"
GTF = "/projects/Annotation/hg19/genecodeV19/gencode.v19.annotation.gtf"
genome_size = "projects/Genome/hg19/hg19.chrom.sizes"



rule all:
    input:
        expand(outdir + "2_rmAdaptor/{sample}_R1.trim.fq.gz", sample=samples),
        expand(outdir + "3_mapping/{sample}.s.bam", sample=samples),
        expand(outdir + "4_rmDup/{sample}.nodup.bam", sample=samples),
        expand(outdir + "4_rmDup/{sample}.nodup.bam.bai", sample=samples),
        expand(outdir + "5_qc_information/{sample}.qc", sample=samples),
        outdir + "5_qc_information/All.qc",
        expand(outdir + "5_bed_bedgraph/{sample}.bedGraph", sample=samples),
        expand(outdir + "5_bed_bedgraph/{sample}.rescale.bw", sample=samples),



#.remove adaptor for pair-end (--pair-filter=any), only keep proper pair reads
rule rm_adaptor:
    input:
        read1= DataDir + "{sample}_R1.fq.gz",
        read2= DataDir + "{sample}_R2.fq.gz"
    log:
        outdir + "2_rmAdaptor/{sample}.log"
    output:
        read1=outdir + "2_rmAdaptor/{sample}_R1.trim.fq.gz",
        read2=outdir + "2_rmAdaptor/{sample}_R2.trim.fq.gz",
    threads: 60
    shell:
        """
        cutadapt -j {threads} -m 18 --quality-cutoff 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --pair-filter=any -o {output.read1} -p {output.read2} {input.read1} {input.read2} 1>{log} 2>&1
        """

rule mapping:
    input:
        fq1 = outdir + "2_rmAdaptor/{sample}_R1.trim.fq.gz",
        fq2 = outdir + "2_rmAdaptor/{sample}_R2.trim.fq.gz",
    output:
        sam = temp(outdir + "3_mapping/{sample}.sam"),
        bam = temp(outdir + "3_mapping/{sample}.bam"),
        Sort = outdir + "3_mapping/{sample}.s.bam",
    params:
        idx = idx,
        unmap = outdir + "3_mapping/{sample}_unmap.fq",
    log: outdir + "3_mapping/{sample}.log",
    threads: 60
    shell:
        """
        bowtie2 -p {threads} -x {params.idx} -1 {input.fq1} -2 {input.fq2} --un-conc {params.unmap} -S {output.sam} --local 1>{log} 2>&1
        samtools view -@ {threads} -bh -o {output.bam} {output.sam}
        samtools sort -@ {threads} -o {output.Sort} {output.bam}
        """

rule picard:
    input:
        bam = outdir + "3_mapping/{sample}.s.bam",
    output:
        bam = outdir + "4_rmDup/{sample}.nodup.bam",
        stat = outdir + "4_rmDup/{sample}.nodup.flagstat",
    log: outdir + "4_rmDup/{sample}.log",
    shell:
        """
        java -jar ~/software/picard.jar MarkDuplicates -I {input.bam} -O {output.bam} \
        --REMOVE_DUPLICATES true -M {log}
    
        samtools flagstat -@ {threads} {output.bam} > {output.stat}
        """

rule makeBai:
    input:
        bam = outdir + "4_rmDup/{sample}.nodup.bam",
    output:
        bai = outdir + "4_rmDup/{sample}.nodup.bam.bai",
    threads: 60
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """

## 统计mapping qc信息
rule qc_count:
    input:
        adaptor = outdir + "2_rmAdaptor/{sample}.log",
        map = outdir + "3_mapping/{sample}.log",
        dup = outdir + "4_rmDup/{sample}.log",
        final = outdir + "4_rmDup/{sample}.nodup.flagstat",
    output:
        outdir + "5_qc_information/{sample}.qc",
    shell:
        """
        name=`basename {input.adaptor} .log`
        totalRead=`grep "Total read pairs processed" {input.adaptor} |cut -d ":" -f 2 |sed 's/ //g'`
        passCutadapt=`grep "Pairs written (passing filters):" {input.adaptor} |cut -d ":" -f 2 |sed 's/ //g'`
        bowtie2_Uniqmap_num=`grep "aligned concordantly exactly 1 time" {input.map} | cut -d " " -f 5`
        bowtie2_Uniqmap_ratio=`grep "aligned concordantly exactly 1 time" {input.map} | cut -d " " -f 6`
        bowtie2_Multimap_num=`grep "aligned concordantly >1 times" {input.map} | cut -d " " -f 5`
        bowtie2_Multimap_ratio=`grep "aligned concordantly >1 times" {input.map} | cut -d " " -f 6`
        PCRdup=`grep "Unknown Library" {input.dup} |awk -F "\\t" '{{print $7"("$9*100"%)"}}'`
        final_tmp=`grep "+ 0 properly paired" {input.final} |cut -d " " -f 1`
        final=`expr $((final_tmp/2))`
        
        echo -e "${{name}}\\t${{totalRead}}\\t${{passCutadapt}}\\t${{bowtie2_Uniqmap_num}}(${{bowtie2_Uniqmap_ratio}})\\t\
        ${{bowtie2_Multimap_num}}(${{bowtie2_Multimap_ratio}})\\t${{PCRdup}}\\t${{final}}" >> {output}
        """

rule cat_qc:
    input:
        expand(outdir + "5_qc_information/{sample}.qc", sample=samples)
    output:
        outdir + "5_qc_information/All.qc",
    shell:
        """
        cat {input} | sed '1i sample\\tRawReads\\trm_Adaptor\\tuniqMap\\tmultiMap\\tPCRdup\\tfinal_R1'> {output}
        """


rule bam2bed2bg:
    input:
        bam = outdir + "4_rmDup/{sample}.nodup.bam",
    output:
        bed = outdir + "5_bed_bedgraph/{sample}.bed",
        bg = outdir + "5_bed_bedgraph/{sample}.bedGraph",
    params:
        chromSize = genome_size,
    shell:
        """
        bedtools bamtobed -i {input.bam} |sort -k1,1 -k2n > {output.bed}
        
        bedtools genomecov -bg -i {output.bed} -g {params.chromSize} > {output.bg}
        """


# python2.7
rule normalize:
    input:
        bed = outdir + "5_bed_bedgraph/{sample}.bed",
        bg = outdir + "5_bed_bedgraph/{sample}.bedGraph",
    output:
        bg = temp(outdir + "5_bed_bedgraph/{sample}.bedGraph.rescale"),
        bgSort = temp(outdir + "5_bed_bedgraph/{sample}.bedGraph.rescale.sort"),
        bw = outdir + "5_bed_bedgraph/{sample}.rescale.bw",
    params:
        chromSize = genome_size,
    shell:
        """
        ~/anaconda3/envs/python2.7/bin/python ~/software/bioutils/scripts/normalize_column.py -c 3 -n 10000000 {input.bg} {input.bed} && \
        
        sort -k1,1 -k2n {output.bg} > {output.bgSort}  && \
        
        bedGraphToBigWig {output.bgSort} {params.chromSize} {output.bw}
        """

