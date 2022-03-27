projectTfDict = {
"BRCA-US-Test": ["FULL"]
}

###  DEFs   ################################################################
# inputs = [expand('{project}_{tfactor}',
    #                  project=key, tfactor=value)
#           for key, value in projectTfDict.items()]

############################################################################


PROJECT = []
TFACTOR = []
for key, value in projectTfDict.items():
    for v in value:
        PROJECT.append(key)
        TFACTOR.append(v)


print(PROJECT)
print(TFACTOR)
############################################################################
rule all:
    input:
        expand("/storage/mathelierarea/processed/petear/analysis/{project}/MethylTable_{project}_{tfactor}.csv.xz", project=PROJECT, tfactor=TFACTOR), \
        expand("/storage/mathelierarea/processed/petear/analysis/{project}/{project}_{tfactor}_universe.bed", project=PROJECT, tfactor=TFACTOR), \
        expand("/storage/mathelierarea/processed/petear/analysis/{project}/Complete_{project}_{tfactor}.pdf", project=PROJECT, tfactor=TFACTOR)

        ############################################################################



# #Rename probes to fit current project
# if TFACTOR == "FULL":
#     rule renameFullProbeslist:
#         input:
#             "FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"
#         output:
#             "{project}_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"
#         shell:
#             "cp *.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed {wildcards.project}_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"

#Remove rows(probes) with more than 50% NAs and impute the rest (Change FROM AMELIA to METHYLIMP when finished):
#metadata with PAM50 and ER.Status is in same location with this name:: sampleinfo_TCGA_RNA_seq_cluster.txt
#metadata with icgc_donor_id etc is at same location with this name::  {project}_sample_Info_260620.tsv
rule FromRawToCis:
    input:
        CpG = "/storage/mathelierarea/processed/petear/SnakemakeInputFiles/{project}_methylation450.sorted.perDonor.tsv",
        coords = "/storage/mathelierarea/processed/petear/SnakemakeInputFiles/{project}_{tfactor}.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed",
        premeta = "/storage/mathelierarea/processed/petear/SnakemakeInputFiles/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/MethylTable_{project}_{tfactor}.csv",
        "/storage/mathelierarea/processed/petear/analysis/{project}/{project}_{tfactor}.pdf",
        "/storage/mathelierarea/processed/petear/analysis/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        "/storage/mathelierarea/processed/petear/analysis/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        "/storage/mathelierarea/processed/petear/analysis/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        "/storage/mathelierarea/processed/petear/analysis/{project}/CTO_{project}_{tfactor}.rds",
        "/storage/mathelierarea/processed/petear/analysis/{project}/binarized_{project}_{tfactor}.csv",
        "/storage/mathelierarea/processed/petear/analysis/{project}/meta_{project}_{tfactor}.csv
    priority:
        100
    script:
        "full.R"

#Python find highest contributing patient and CpG to topics and make violinplots with seaborn.
rule topPatientCpG:
    input:
        methTab = "/storage/mathelierarea/processed/petear/analysis/{project}/MethylTable_{project}_{tfactor}.csv",
        regScrNorm = "/storage/mathelierarea/processed/petear/analysis/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        regScrUnrm = "/storage/mathelierarea/processed/petear/analysis/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        topicAssig = "/storage/mathelierarea/processed/petear/analysis/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        binFrame = "/storage/mathelierarea/processed/petear/analysis/{project}/binarized_{project}_{tfactor}.csv"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/Violin_{project}_{tfactor}.pdf"
        "/storage/mathelierarea/processed/petear/analysis/{project}/ProbeTopicScore_{project}_{tfactor}.csv"
    script:
        "seabornViolin.py"


#Make heatmap in ggplot2 from top patients and CpGs found in 'rule seabornViolin'
rule probeHeatmap:
    input:
        "/storage/mathelierarea/processed/petear/analysis/{project}/ProbeTopicScore_{project}_{tfactor}.csv"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/PrbTpcScrHeatmap_{project}_{tfactor}.pdf"
    script:
        "Heatmap.R"

#CSVs can be big so this rule compresses the CSV to XZ to save space.
rule MethTabl_CSVtoXZ:
    input:
        CSV = "/storage/mathelierarea/processed/petear/analysis/{project}/MethylTable_{project}_{tfactor}.csv"
    output:
        XZ = "/storage/mathelierarea/processed/petear/analysis/{project}/MethylTable_{project}_{tfactor}.csv.xz"
    shell:
        "xz --threads=0 -9v {input.CSV}"


#Make outputfolder and add variable number of bedfiles there:
checkpoint get_bed:
    input:
        "/storage/mathelierarea/processed/petear/analysis/{project}/CTO_{project}_{tfactor}.rds"
    output:
        directory("/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}")
    priority:
        99
    shell:
        "Rscript /storage/mathelierarea/processed/petear/Snakemake2/CisBed.R /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/CTO_{wildcards.project}_{wildcards.tfactor}.rds /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}"


#Liftover bed from HG19 to HG38 and then delete the HG19 (To easier use snakemake with unibind)
rule liftover:
    input:
        HG19Bed = "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}.bed"
    output:
        HG38Bed = "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    priority:
        98
    shell:
        "liftOver {input.HG19Bed} /storage/mathelierarea/processed/petear/hg19ToHg38.over.chain.gz {output.HG38Bed} unMapped"





#make function to continue to work on sepearate bedfiles
def cont_work_bed(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names


#Create universe (Background)
rule universe:
    input:
        cont_work_bed
    output:
        Universe = "/storage/mathelierarea/processed/petear/analysis/{project}/{project}_{tfactor}_universe.bed"
    shell:
        "cat /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/*_HG38.bed | bedtools sort | bedtools merge > {output}"


#Run unibind enrichment analysis on each bed file from checkpoint get_bed. (Should this also be a checkpoint?) .
rule unibind:
    input:
        topicBed = "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
        Universe = "/storage/mathelierarea/processed/petear/analysis/{project}/{project}_{tfactor}_universe.bed"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf"

    shell:
        "bash /storage/mathelierarea/processed/petear/UnibindFolder/UnibindMaster/bin/UniBind_enrich.sh oneSetBg /storage/mathelierarea/processed/petear/UnibindFolder/LolaUpdated/UniBind_LOLA.RDS {input.topicBed} {input.Universe} /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/UB_output_Topic{wildcards.bednumber}"


#make function to continue to work on sepearate bedfiles
def unibindFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names

rule mergeUBPDF:
    input:
        unibindFunc
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"

rule mergeCisTPDF:
    input:
        cistopicPDF = "/storage/mathelierarea/processed/petear/analysis/{project}/{project}_{tfactor}.pdf",
        UBpdf = "/storage/mathelierarea/processed/petear/analysis/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/merged_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input.cistopicPDF} {input.UBpdf} {output}"

rule testFunc:
    input:
        unibindFunc
    output:
        "test_{project}_{tfactor}.txt"
    shell:
        "echo 'done' > {output}"

# #&& rm Topic_{wildcards.bednumber}.bed

    ########################################
rule rGREAT:
    input:
        HG38Bed = "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    output:
        GreatPDF = "/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf"
    shell:
        "Rscript /storage/mathelierarea/processed/petear/Snakemake2/GREAT.R /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/Topic_{wildcards.bednumber}_HG38.bed /storage/mathelierarea/processed/petear/analysis/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/{wildcards.bednumber}_HG38.bed_GREAT.pdf"

def rGREATFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/storage/mathelierarea/processed/petear/analysis/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names


rule MergeGREATPDFs:
    input:
        rGREATFunc
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"


rule completePDF:
    input:
        cisTPDF =   "/storage/mathelierarea/processed/petear/analysis/{project}/merged_{project}_{tfactor}.pdf",
        GreatPDF = "/storage/mathelierarea/processed/petear/analysis/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    output:
        "/storage/mathelierarea/processed/petear/analysis/{project}/Complete_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input.cisTPDF} {input.GreatPDF} {output}"
