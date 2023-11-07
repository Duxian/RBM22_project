##################################################
# GRO-seq data analysis pipline 
##################################################

configfile: "config5.yaml"

SAMPLES = config["SAMPLES"]
DataDir = config["DataDir"]
outdir = config["outdir"]

genome_idx = "/projects/Genome/star_idx/hg19"
GTF = "/projects/Annotation/hg19/genecodeV19/gencode.v19.annotation.gtf"
genome_size = "/projects/Genome/hg19/hg19.chrom.sizes"
genesModel = "/projects/Annotation/hg19/genecodeV19/annoFile/Genes.bed"



rule all:
    input:
        expand(outdir + "1_fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand(outdir + "2_cutadapt/{sample}_R1.trim.fq.gz", sample=SAMPLES),
        expand(outdir + "3_fastqc2/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand(outdir + "4a_rm_rRNA/{sample}.sam", sample=SAMPLES),
        expand(outdir + "4b_star_hg19/{sample}Aligned.out.bam", sample=SAMPLES),
        expand(outdir + "5_BamSort/uniqMap/{sample}.uniq.bed", sample=SAMPLES),
        expand(outdir + "5_BamSort/uniqMap/{sample}.strand.log", sample=SAMPLES),
        expand(outdir + "6_Track/bw/{sample}.Forward.bw", sample=SAMPLES),
        expand(outdir + "7_qc_information/{sample}.qc", sample=SAMPLES),
        expand(outdir + "7_qc_information/All.qc", sample=SAMPLES),


#1.fastqc round1 to check quality and adaptor
rule fastqc_round1:
    input:
        read1 = DataDir + "{sample}_R1.fq.gz",
        read2 = DataDir + "{sample}_R2.fq.gz",
    log:
        log1 = outdir + "1_fastqc/{sample}_R1.log",
        log2 = outdir + "1_fastqc/{sample}_R2.log",
    output:
        html1 = outdir + "1_fastqc/{sample}_R1_fastqc.html",
        html2 = outdir + "1_fastqc/{sample}_R2_fastqc.html",
    params:
        outdir + "1_fastqc/"
    threads: 60
    shell:
        """
        fastqc {input.read1} -o {params} -t {threads} > {log.log1} 2>&1
        fastqc {input.read2} -o {params} -t {threads} > {log.log2} 2>&1
        """

rule rm_adaptor:
    input:
        read1= DataDir + "{sample}_R1.fq.gz",
        read2= DataDir + "{sample}_R2.fq.gz"
    log:
        outdir + "2_cutadapt/{sample}.log"
    output:
        read1=outdir + "2_cutadapt/{sample}_R1.trim.fq.gz",
        read2=outdir + "2_cutadapt/{sample}_R2.trim.fq.gz",
    threads: 60
    shell:
        """
        cutadapt -j {threads} -m 40 --quality-cutoff 20 -a AGATCGGAAGAGCACACG -A AGATCGGAAGAGCGTCGTG \
        --pair-filter=any -o {output.read1} -p {output.read2} {input.read1} {input.read2} 1>{log} 2>&1
        """


#3.fastqc round2 to check remove adaptor results
rule fatqc_round2:
    input:
        read1 = outdir + "2_cutadapt/{sample}_R1.trim.fq.gz",
        read2 = outdir + "2_cutadapt/{sample}_R2.trim.fq.gz",
    log:
        outdir + "3_fastqc2/{sample}.log"
    output:
        html1 = outdir + "3_fastqc2/{sample}_R1_fastqc.html",
        html2 = outdir + "3_fastqc2/{sample}_R2_fastqc.html"
    params:
        outdir + "3_fastqc2/"
    threads: 60
    shell:
        """
        fastqc {input.read1} -o {params} -t {threads} > {log} 2>&1
        fastqc {input.read2} -o {params} -t {threads} > {log} 2>&1
        """

rule STAR_map_hg19:
    input:
        read1 = outdir + "2_cutadapt/{sample}_R1.trim.fq.gz",
    params:
        prefix = outdir + "4b_star_hg19/{sample}",
        GTF = GTF,
        genome_idx = genome_idx,
    log:
        outdir + "4b_star_hg19/{sample}.log"
    threads: 60
    output:
        bam = outdir + "4b_star_hg19/{sample}Aligned.out.bam",
        log = outdir + "4b_star_hg19/{sample}Log.final.out",
    shell:
        """
	    STAR --runMode alignReads --runThreadN {threads} \
        --genomeDir {params.genome_idx} --sjdbGTFfile {params.GTF} \
        --genomeLoad NoSharedMemory \
	    --sjdbOverhang 149 \
        --readFilesIn {input.read1}  --readFilesCommand gunzip -c \
        --alignEndsType Local \
	    --outFilterMismatchNoverLmax 0.03  \
        --outFilterMatchNminOverLread 0.5 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.5 \
        --outSAMunmapped None --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outFileNamePrefix {params.prefix}  \
        --outSAMtype BAM Unsorted --outFilterType BySJout --outSAMattributes All \
        --outReadsUnmapped Fastx 1>{log} 2>&1        
	    """

rule picard:
    input:
        bam = outdir + "4b_star_hg19/{sample}Aligned.out.bam",
    output:
        sort =  temp(outdir + "4b_star_hg19/{sample}Aligned.out.s.bam"),
        bam = outdir + "4b_star_hg19/{sample}.nodup.bam",
    log: outdir + "4b_star_hg19/{sample}_markdup_metrics.txt",
    threads: 60
    shell:
        """
        samtools sort -o {output.sort} -@ {threads} {input.bam}
        
        java -jar ~/software/picard.jar MarkDuplicates -I {output.sort} -O {output.bam} --REMOVE_DUPLICATES true -M {log}
        """

rule uniqBam:
    input:
        bam = outdir + "4b_star_hg19/{sample}.nodup.bam",
    output:
        bam = temp(outdir + "5_BamSort/uniqMap/{sample}.uniq.bam"),
        sortBam = outdir + "5_BamSort/uniqMap/{sample}.uniq.s.bam",
        bed = outdir + "5_BamSort/uniqMap/{sample}.uniq.bed",
        qc = outdir + "5_BamSort/uniqMap/{sample}.final.qc",
    threads: 60
    shell:
        """
        samtools view -@ {threads} -q 255 -b -h {input.bam} -o {output.bam}
        samtools sort -o {output.sortBam} -@ {threads} {output.bam}
        samtools index -@ {threads} {output.sortBam} 
        bedtools bamtobed -i {output.sortBam} -bed12 -split > {output.bed}
        
        samtools flagstat -@ {threads} {output.sortBam} > {output.qc}
        """

# guess strand
rule judge_strand:
    input:
        bamsort = outdir + "5_BamSort/uniqMap/{sample}.uniq.s.bam",
    output:
        strandLog = outdir + "5_BamSort/uniqMap/{sample}.strand.log",
    params:
        genesModel = genesModel
    shell:
        """
        infer_experiment.py -r {params.genesModel} -i {input.bamsort} > {output.strandLog}
        """

## 统计mapping qc信息
rule qc_count:
    input:
        adaptor = outdir + "2_cutadapt/{sample}.log",
        map = outdir + "4b_star_hg19/{sample}Log.final.out",
        dup = outdir + "4b_star_hg19/{sample}_markdup_metrics.txt",
        final = outdir + "5_BamSort/uniqMap/{sample}.final.qc",
    output:
        outdir + "7_qc_information/{sample}.qc",
    shell:
        """
        name=`basename {input.adaptor} .log`
        totalRead=`grep "Total read pairs processed" {input.adaptor} |cut -d ":" -f 2 |sed 's/ //g'`
        passCutadapt=`grep "Pairs written (passing filters):" {input.adaptor} |cut -d ":" -f 2 |sed 's/ //g'`
        star_Uniqmap_num=`grep "Uniquely mapped reads number" {input.map} | cut -f 2`
        star_Uniqmap_ratio=`grep "Uniquely mapped reads %"  {input.map} | cut -f 2`
        star_Multimap_num=`grep "Number of reads mapped to multiple loci" {input.map} | cut -f 2`
        star_Multimap_ratio=`grep "% of reads mapped to multiple loci" {input.map} | cut -f 2`
        PCRdup=`grep "Unknown Library" {input.dup} |awk -F "\\t" '{{print $6"("$9*100"%)"}}'`
	    final=`grep "mapped (100.00% : N/A)" {input.final} |head -1 |cut -d " " -f 1`
        
        echo -e "${{name}}\\t${{totalRead}}\\t${{passCutadapt}}\\t${{star_Uniqmap_num}}(${{star_Uniqmap_ratio}})\\t${{star_Multimap_num}}(${{star_Multimap_ratio}})\\t${{PCRdup}}\\t${{final}}" >> {output}
        """

rule cat_qc:
    input:
        expand(outdir + "7_qc_information/{sample}.qc", sample=SAMPLES)
    output:
        outdir + "7_qc_information/All.qc",
    shell:
        """
        cat {input} | sed '1i sample\\tRawReads\\trm_Adaptor\\tuniqMap\\tmultiMap\\tPCRdup\\tfinal'> {output}
        """


# bam to bigiwig and not normalize
rule bam2wig:
    input:
        bamsort = outdir + "5_BamSort/uniqMap/{sample}.uniq.s.bam",
    output:
        bigwig = outdir + "6_Track/bw/{sample}.Forward.bw",
    params:
        genome_size = genome_size,
        prefix = outdir + "6_Track/bw/{sample}",
    threads: 60
    shell:
        """
        bam2wig.py -i {input.bamsort} -s {params.genome_size} -o {params.prefix} --strand='+-,-+' -t 1000000000
        """

