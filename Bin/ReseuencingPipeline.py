#!/usr/bin/env python2

import sys
import re
import os
import argparse
import ConfigParser
from glob import glob

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
FASTP = config.get('SOFTWARE', 'fastp')
BWA = config.get('SOFTWARE', 'bwa')
SAMTOOLS = config.get('SOFTWARE', 'samtools')
GATK = config.get('SOFTWARE', 'gatk4')
BEDTOOLS = config.get('SOFTWARE', 'bedtools')
PYTHON = config.get('SOFTWARE', 'python')
SNAKEMAKE = config.get('SOFTWARE', 'snakemake')

### SCRIPT

#### DATABASE

class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\s+', line_in)
        self.Sample = list_split[0]
        self.Group = list_split[1]
        self.Name = '_'.join(list_split[:2])
        if ',' in list_split[2]:
            self.paired = True
            list_fq = re.split(',', list_split[2])
            self.fq1 = list_fq[0]
            self.fq2 = list_fq[1]
        else:
            self.paired = False
            self.fq = list_split[2]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')


def Argparse():
    parser = argparse.ArgumentParser(description="PAS pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the abs output path', required=True)
    parser.add_argument('-a1', help='the read1 adapter', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-a2', help='the read2 adapter', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-p', help='the name of spe', default='hg19')
    parser.add_argument('-t', help='the core used', default=10)
    parser.add_argument('-sge', help='use sge',  action='store_true')
    parser.add_argument('-lsf', help='use lsf',  action='store_true')
    parser.add_argument('-r', help='run now',  action='store_true')
    argv=vars(parser.parse_args())
    return argv

def HandleRawdata(argv):
    outpath = argv['o']
    RawData = os.path.join(outpath, '0.RawData')
    list_ob = []
    if not os.path.exists(RawData):
        os.mkdir(RawData)
    else:
        os.system('rm -rf '+RawData)
        os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                if ob.paired:
                    Paired = True
                    if ob.fq1.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq1, ob.Name+'_1'+lb))
                    os.system('ln -s {} {}'.format(ob.fq2, ob.Name+'_2'+lb))
                else:
                    Paired = False
                    if ob.fq.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq, ob.Name+lb))
    os.chdir(outpath)
    return list_ob, Paired, lb

def WriteSnake(argv, list_ob, Paired, lb):
    DB = config.get('DATABASE', argv['p'])
    outpath = argv['o']
    label = os.path.basename(glob(DB+'/*.fa')[0]).rstrip('.fa')
    snakefile = open(os.path.join(argv['o'], 'snakefile.txt'), 'w')
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))
    snakefile.write('adapter1 = "{}"\nadapter2 = "{}"\n'.format(argv['a1'], argv['a2']))
    snakefile.write('FASTA = "{}"\n'.format(DB+'/'+label+'.fa'))
    snakefile.write('REGION = "{}"\n'.format(DB+'/'+label+'.fa.bed'))

    ###all
    All = Snake('All')
    All.UpdateInput('"4.Mutation/All.raw.filter.vcf"')
    All.WriteStr(snakefile)

    ###QC
    QC = Snake('QC')
    if Paired:
        QC.UpdateInput('A = "'+outpath+'/0.RawData/{sample}_1'+lb+'", B = "'+outpath+'/0.RawData/{sample}_2'+lb+'"')
        QC.UpdateOutput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input.A} -o {output.A} -I {input.B} -O {output.B} --adapter_sequence {adapter1} --adapter_sequence_r2 {adapter2} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    else:
        QC.UpdateInput('"'+outpath+'/0.RawData/{sample}'+lb+'"')
        QC.UpdateOutput('"'+outpath+'/1.QC/{sample}.clean.fastq.gz"')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input} -o {output}  --adapter_sequence {adapter1} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html"')
    QC.WriteStr(snakefile)

    ### Align
    Align = Snake('Align')
    if Paired:
        Align.UpdateInput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        Align.UpdateOutput('"2.Align/{sample}.bam"')
        Align.UpdateThreads('5')
        Align.UpdateLog('e = "logs/{sample}.align.e", o = "logs/{sample}.align.o"')
        Align.UpdateParams('rg="@RG\\\\tID:{sample}\\\\tPL:illumina\\\\tSM:{sample}"')
        Align.UpdateShell(r'"'+BWA+r" mem -t 5 -M -R '{params.rg}' {FASTA} {input.A} {input.B} |"+SAMTOOLS+r' view -Sb -@ 4 -T {FASTA} -o {output}"')
    else:
        Align.UpdateInput('"'+outpath+'/1.QC/{sample}.clean.fastq.gz"')
        Align.UpdateOutput('"2.Align/{sample}.bam"')
        Align.UpdateThreads('5')
        Align.UpdateLog('e = "logs/{sample}.align.e", o = "logs/{sample}.align.o"')
        Align.UpdateParams('rg="@RG\\\\tID:{sample}\\\\tPL:illumina\\\\tSM:{sample}"')
        Align.UpdateShell(r'"'+BWA+r" mem -t 5 -M -R '{params.rg}' {FASTA} {input} |"+SAMTOOLS+r' view -Sb -@ 4 -T {FASTA} -o {output}"')
    Align.WriteStr(snakefile)

    ### Sort
    Sort = Snake('Sort')
    Sort.UpdateInput('"2.Align/{sample}.bam"')
    Sort.UpdateOutput('"2.Align/{sample}.sort.bam"')
    Sort.UpdateLog('e = "logs/{sample}.sort.e", o = "logs/{sample}.sort.o"')
    Sort.UpdateShell(r'"'+GATK+r' SortSam --INPUT {input} --OUTPUT {output} --SORT_ORDER coordinate"')
    Sort.WriteStr(snakefile)

    ### MarkDup
    MarkDup = Snake('MarkDup')
    MarkDup.UpdateInput('"2.Align/{sample}.sort.bam"')
    MarkDup.UpdateOutput('"2.Align/{sample}.sort.dedup.bam"')
    MarkDup.UpdateLog('e = "logs/{sample}.redup.e", o = "logs/{sample}.redup.o"')
    MarkDup.UpdateShell(r'"'+GATK+r' MarkDuplicates --INPUT {input} --OUTPUT {output} --METRICS_FILE 2.Align/{wildcards.sample}.metrics --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000"')
    MarkDup.WriteStr(snakefile)

    ### BamIndex
    BamIndex = Snake('BamIndex')
    BamIndex.UpdateInput('"2.Align/{sample}.sort.dedup.bam"')
    BamIndex.UpdateOutput('"2.Align/{sample}.sort.dedup.bam.bai"')
    BamIndex.UpdateLog('e = "logs/{sample}.BamIndex.e", o = "logs/{sample}.BamIndex.o"')
    BamIndex.UpdateShell(r'"'+SAMTOOLS+r' index {input}"')
    BamIndex.WriteStr(snakefile)
    '''

    ### BaseRecalibrator
    BaseRecalibrator = Snake('BaseRecalibrator')
    BaseRecalibrator.UpdateInput('"2.Align/{sample}.sort.dedup.bam.bai"')
    BaseRecalibrator.UpdateOutput('"3.Call/{sample}.recalibration.report"')
    BaseRecalibrator.UpdateLog('e = "logs/{sample}.BaseRecalibrator.e", o = "logs/{sample}.BaseRecalibrator.o"')
    BaseRecalibrator.UpdateShell(r'"'+GATK+' BaseRecalibrator --input 2.Align/{wildcards.sample}.sort.dedup.bam --output {output} --reference {FASTA} --intervals {REGION}"')
    BaseRecalibrator.WriteStr(snakefile)

    ### GatherBQSRReports
    GatherBQSRReports = Snake('GatherBQSRReports')
    GatherBQSRReports.UpdateInput('expand("3.Call/{sample}.recalibration.report", sample=Samples)')
    insh = '-I '+' -I '.join(["3.Call/"+i.Name+".recalibration.report" for i in list_ob])
    GatherBQSRReports.UpdateOutput('"3.Call/all.recalibration.report"')
    GatherBQSRReports.UpdateLog('e = "logs/GatherBQSRReports.e", o = "logs/GatherBQSRReports.o"')
    GatherBQSRReports.UpdateShell(r'"'+GATK+' GatherBQSRReports '+insh+' -O {output}"')
    GatherBQSRReports.WriteStr(snakefile)

    ### ApplyBQSR
    ApplyBQSR = Snake('ApplyBQSR')
    ApplyBQSR.UpdateInput('report="3.Call/all.recalibration.report",bam="2.Align/{sample}.sort.dedup.bam"')
    ApplyBQSR.UpdateOutput('"3.Call/{sample}.sort.dedup.BQSR.bam"')
    ApplyBQSR.UpdateLog('e = "logs/{sample}.ApplyBQSR.e", o = "logs/{sample}.ApplyBQSR.o"')
    ApplyBQSR.UpdateShell(r'"'+GATK+' ApplyBQSR -R {FASTA} -I {input.bam} -O {output} -bqsr {input.report} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --use-original-qualities -L {REGION}"')
    ApplyBQSR.WriteStr(snakefile)
    '''
    ### HaplotypeCaller
    HaplotypeCaller = Snake('HaplotypeCaller')
    HaplotypeCaller.UpdateInput('A = "2.Align/{sample}.sort.dedup.bam", B = "2.Align/{sample}.sort.dedup.bam.bai"')
    HaplotypeCaller.UpdateOutput('"3.Call/{sample}.gvcf"')
    HaplotypeCaller.UpdateLog('e = "logs/{sample}.HaplotypeCaller.e", o = "logs/{sample}.HaplotypeCaller.o"')
    HaplotypeCaller.UpdateShell(r'"'+GATK+' HaplotypeCaller -R {FASTA} --emit-ref-confidence GVCF -I {input.A} -L {REGION} -O {output}"')
    HaplotypeCaller.WriteStr(snakefile)

    ### CombineGVCFs
    CombineGVCFs = Snake("CombineGVCFs")
    cmsh = '-V '+' -V '.join(["3.Call/"+i.Name+".gvcf" for i in list_ob])
    CombineGVCFs.UpdateInput('expand("3.Call/{sample}.gvcf", sample=Samples)')
    CombineGVCFs.UpdateOutput('"3.Call/All.gvcf"')
    CombineGVCFs.UpdateLog('e = "logs/CombineGVCFs.e", o = "logs/CombineGVCFs.o"')
    CombineGVCFs.UpdateShell(r'"'+GATK+' CombineGVCFs -R {FASTA} '+cmsh+' -O {output}"')
    CombineGVCFs.WriteStr(snakefile)

    ### GenotypeGVCFs
    GenotypeGVCFs = Snake("GenotypeGVCFs")
    GenotypeGVCFs.UpdateInput('"3.Call/All.gvcf"')
    GenotypeGVCFs.UpdateOutput('"3.Call/All.raw.gvcf"')
    GenotypeGVCFs.UpdateLog('e = "logs/GenotypeGVCFs.e", o = "logs/GenotypeGVCFs.o"')
    GenotypeGVCFs.UpdateShell(r'"'+GATK+' GenotypeGVCFs -R {FASTA} -V {input} -O {output} -L {REGION}"')
    GenotypeGVCFs.WriteStr(snakefile)

    ### SelectSNP
    SelectSNP = Snake("SelectSNP")
    SelectSNP.UpdateInput('"3.Call/All.raw.gvcf"')
    SelectSNP.UpdateOutput('"4.Mutation/All.raw.snp.vcf"')
    SelectSNP.UpdateLog('e = "logs/SelectSNP.e", o = "logs/SelectSNP.o"')
    SelectSNP.UpdateShell(r'"'+GATK+' SelectVariants -select-type SNP -V {input} -O {output}"')
    SelectSNP.WriteStr(snakefile)

    ### SelectIndel
    SelectIndel = Snake("SelectIndel")
    SelectIndel.UpdateInput('"3.Call/All.raw.gvcf"')
    SelectIndel.UpdateOutput('"4.Mutation/All.raw.indel.vcf"')
    SelectIndel.UpdateLog('e = "logs/SelectIndel.e", o = "logs/SelectIndel.o"')
    SelectIndel.UpdateShell(r'"'+GATK+' SelectVariants -select-type INDEL -V {input} -O {output}"')
    SelectIndel.WriteStr(snakefile)

    ### FilterSNP
    FilterSNP = Snake("FilterSNP")
    FilterSNP.UpdateInput('"4.Mutation/All.raw.snp.vcf"')
    FilterSNP.UpdateOutput('"4.Mutation/All.raw.snp.filter.vcf"')
    FilterSNP.UpdateLog('e = "logs/FilterSNP.e", o = "logs/FilterSNP.o"')
    FilterSNP.UpdateShell('\''+GATK+' VariantFiltration -V {input} --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O {output}'+'\'')
    FilterSNP.WriteStr(snakefile)

    ### FilterIndel
    FilterIndel = Snake("FilterIndel")
    FilterIndel.UpdateInput('"4.Mutation/All.raw.indel.vcf"')
    FilterIndel.UpdateOutput('"4.Mutation/All.raw.indel.filter.vcf"')
    FilterIndel.UpdateLog('e = "logs/FilterIndel.e", o = "logs/FilterIndel.o"')
    FilterIndel.UpdateShell('\''+GATK+' VariantFiltration -V {input} --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O {output}'+'\'')
    FilterIndel.WriteStr(snakefile)

    ### MergeVcfs
    MergeVcfs = Snake("MergeVcfs")
    MergeVcfs.UpdateInput('indel="4.Mutation/All.raw.indel.filter.vcf",snp="4.Mutation/All.raw.snp.filter.vcf"')
    MergeVcfs.UpdateOutput('"4.Mutation/All.raw.filter.vcf"')
    MergeVcfs.UpdateLog('e = "logs/MergeVcfs.e", o = "logs/MergeVcfs.o"')
    MergeVcfs.UpdateShell(r'"'+GATK+' MergeVcfs -I {input.snp} -I {input.indel} -O {output}"')
    MergeVcfs.WriteStr(snakefile)

def RunShell(argv):
    out = open(os.path.join(argv['o'], 'runsnake.sh'), 'w')
    if argv['lsf']:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--cluster', "'bsub -q normal -n {threads} -o %J.o -e %J.e'",\
                                      '--printshellcmds', '--snakefile', 'snakefile.txt']))
    else:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--printshellcmds', '--snakefile', 'snakefile.txt']))
    out.close()
    if argv['r']:
        os.system('nohup sh runsnake.sh& > snake.log')

def main():
    argv = Argparse()
    list_ob, Paired, lb = HandleRawdata(argv)
    WriteSnake(argv, list_ob, Paired, lb)
    RunShell(argv)


if __name__ == '__main__':
    main()
