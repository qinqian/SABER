import glob
import os.path
import re
configfile: "config.yaml"

adapters = ['XNNNNNNNNNNACTCAGATCTCGAGCTCAAGCTTCG', 
            'GAGACAAATGGCAGCCCG']

barcode = 'NNNNNNNNNN'


def INPUTFOLDER(wildcards):
    output = config['input_prefix']
    if wildcards.project == "testing":
        output="testing_"+output
    elif wildcards.project == "validation":
        output="validation_"+output
    elif wildcards.project == "demo":
        output="demo"
    return output

def INPUTSAMPLES(wildcards):
    output = INPUTFOLDER(wildcards)
    print(output)
    fq_files=glob.glob(os.path.join(output, "*.fastq.gz"))
    print(fq_files)
    r = re.compile("(.*)(_S[0-9].*)?(_L00[1-9])?_R[12](_001)?.fastq.gz")
    fq_files= list(filter(r.match,fq_files))
    SAMPLES=list(set([r.match(x).group(1)[len(output)+1::] for x in fq_files]))
    return SAMPLES

def remapSampleToFastq(wildcards):
    # Miseq files end in _001.fastq.gz ; NWCH Hiseq files do not
    #r1list = glob.glob(os.path.join(config['input_prefix'], wildcards.sample + "_*_R1.fastq.gz")) # usu ends in _001.fastq.gz; not NWCH
    #r2list = glob.glob(os.path.join(config['input_prefix'], wildcards.sample + "_*_R2.fastq.gz"))
    r1_glob_pattern = os.path.join(INPUTFOLDER(wildcards), wildcards.sample + "_*R1")
    r2_glob_pattern = os.path.join(INPUTFOLDER(wildcards), wildcards.sample + "_*R2")
    r1_glob_pattern += "*.fastq.gz"
    r2_glob_pattern += "*.fastq.gz"

    r1list = glob.glob(r1_glob_pattern)
    r2list = glob.glob(r2_glob_pattern)

    r1file = r1list[0]
    r2file = r2list[0]
    return [r1file, r2file]

def umi_trim_switch(wildcards):
    if(config["use_umi"]):
        return wildcards.project+"_output/fastq_umi/"+wildcards.sample+".fastq"
    else:
        return wildcards.project+"_output/cutadapt3p/"+wildcards.sample+".fastq" # trim too much
        #return wildcards.project+"_output/fastq_trimmed/"+wildcards.sample+".trimmed.fastq"

def umi_bam_switch(wildcards):
    if(config["use_umi"]):
        return expand(wildcards.project+"_output/dedup/{sample}.bam.bai",sample=INPUTSAMPLES(wildcards))
    else:
        return expand(wildcards.project+"_output/bam/{sample}.bam.bai",sample=INPUTSAMPLES(wildcards))

def multiqc_input(wildcards):
    return expand("{project}_output/fastqc/{sample}_fastqc.html",sample=INPUTSAMPLES(wildcards),project=wildcards.project)

rule normal_run:
    input: config["output_prefix"]+"_output/SABER/analysis.done",config["output_prefix"]+"_output/multiqc/multiqc_report.html"

rule testing:
    input: "testing_output/SABER/analysis.done","testing_output/multiqc/multiqc_report.html"

rule validation:
    input: "validation_output/SABER/analysis.done","validation_output/multiqc/multiqc_report.html"

rule demo:
    input: "demo_output/SABER/analysis.done","demo_output/multiqc/multiqc_report.html"

rule decompress:
    input:"%s/{sample}.fastq.gz" % config['input_prefix']
    output:"%s/{sample}.fastq" % config['input_prefix']
    # TODO: should I assume gzip is installed?
    shell:"gzip -d -c {input} > {output}"

rule merge_read_pairs:
    input:remapSampleToFastq
    output:"{project}_output/fastq_merged/{sample}.assembled.fastq"
    log:"{project}_logs/pear/{sample}.PEARreport.txt"
    conda:"envs/tools-env.yaml"
    threads:1
    shell:"pear -j {threads} -f {input[0]} -r {input[1]} -o {wildcards.project}_output/fastq_merged/{wildcards.sample} > {log}"

rule trimmomatic:
    input:"{project}_output/fastq_merged/{sample}.assembled.fastq"
    output:"{project}_output/fastq_trimmed/{sample}.trimmed.fastq"
    conda:"envs/tools-env.yaml"
    shell:"trimmomatic SE {input} {output} SLIDINGWINDOW:4:15 MINLEN:100"

#run cutadapt to keep only fastqs that have the full flanking primer sequences
#5 prime adapter = V6/7F:  TCGAGCTCAAGCTTCGG
rule cutadapt5p:
    input:"{project}_output/fastq_trimmed/{sample}.trimmed.fastq"
    output:"{project}_output/cutadapt5p/{sample}.fastq"
    params:
        #adapter= "XCGCAGAGAGGCTCCGTG" if config["use_umi"] else "XTCGAGCTCAAGCTTCGG"
        adapter= adapters[0]
        #adapter= "XCTAAATGGCTGTGAGAGAGCTCAG" if config["use_umi"] else "XACTCAGATCTCGAGCTCAAG"
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell:"cutadapt -g {params.adapter} --discard-untrimmed -e 0.01 --action=none -o {output} {input}"

#3 prime adapter = V6/7R:  GACCTCGAGACAAATGGCAG (reverse complement of the primer sequence 5'-3')
rule cutadapt3p:
    input:"{project}_output/cutadapt5p/{sample}.fastq"
    output:"{project}_output/cutadapt3p/{sample}.fastq"
    conda:"envs/tools-env.yaml"
    #changed to allow full adapter sequence anywhere to accomodate UMI.  If this doesn't work change X to ^.
    shell: "cutadapt -a %s$ --discard-untrimmed -e 0.01 --action=none -o {output} {input}" % adapters[1] ## working one 


# run fastqc on trimmed and demultiplexed input files
rule fastqc:
    input:"{project}_output/cutadapt3p/{sample}.fastq"
    output:"{project}_output/fastqc/{sample}_fastqc.html","{project}_output/fastqc/{sample}_fastqc.zip"
    conda:"envs/tools-env.yaml"
    shell:"fastqc -o {wildcards.project}_output/fastqc/ {input}"

# run multiqc to summarize the qc files
rule multiqc:
    input: multiqc_input
    output:"{project}_output/multiqc/multiqc_report.html"
    conda:"envs/tools-env.yaml"
    shell:"multiqc -d {wildcards.project}_output/fastqc/ -o {wildcards.project}_output/multiqc/"

# extract UMI tags (optional)
rule umi_extract:
    input:"{project}_output/cutadapt3p/{sample}.fastq"
    output:"{project}_output/fastq_umi/{sample}.fastq"
    log:"{project}_logs/umi_extract/{sample}.log"
    conda:"envs/tools-env.yaml"
    #shell:"umi_tools extract --stdin={input} --bc-pattern=XXXXXXXXXXXXXXXXXNNNNNNNNNN --log={log} --stdout={output}"
    shell:"umi_tools extract --stdin={input} --bc-pattern=%s --log={log} --stdout={output}" % barcode

rule needleall:
    input:umi_trim_switch
    output:"{project}_output/needleall/{sample}.sam"
    log:"{project}_logs/needleall/{sample}.error"
    conda:"envs/tools-env.yaml"
    shell:"needleall -aformat3 sam -gapextend 0.25 -gapopen 10.0 -awidth3=5000 -asequence references/hspGESTALT_amplicon_corrected.fa -bsequence {input} -outfile {output} -errfile {log}"

# fix sam file header and select only reads matching at position 1
rule altersam:
    input:
        sam="{project}_output/needleall/{sample}.sam",
        header="references/hspGESTALT_header.tsv"
    output:"{project}_output/bam/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"cat {input.header} <(grep -v ^@ {input.sam} | awk '{{if ($4==1){{print}}}}') | samtools view -hb | samtools sort> {output}"

rule index:
    input:"{fn}.bam"
    output:"{fn}.bam.bai"
    conda:"envs/tools-env.yaml"
    shell:"samtools index {input}"

rule dedup:
    input:
        bam="{project}_output/bam/{sample}.bam",
        bai="{project}_output/bam/{sample}.bam.bai"
    output:"{project}_output/dedup/{sample}.bam"
    conda:"envs/tools-env.yaml"
    shell:"umi_tools dedup --method=unique -I {input.bam} --output-stats= output/dedup/{wildcards.sample} -S {output}"

rule analysis:
    input:umi_bam_switch
    output:touch("{project}_output/SABER/analysis.done")
    conda:"envs/r-env.yaml"
    params:
        bams= lambda wildcards: wildcards.project+"_output/dedup" if config["use_umi"] else wildcards.project+"_output/bam",
        out= lambda wildcards: wildcards.project+"_output/SABER",
        conf= "config.yaml"
    shell:"""
    Rscript scripts/installrpy.R || :
    ##Rscript scripts/analysis.R {params.conf} {params.bams} {params.out}
    Rscript scripts/analysis_customize.R {params.conf} {params.bams} {params.out}
    """
