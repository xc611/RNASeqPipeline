from snakemake.utils import R

configfile: "config.json"

# output and input directory

workdir: config["MYFOLDER"]

INPUTDIR= config["INPUTDIR"]

## sample names and contrast information

SAMPLES= config["MYSAMPLES"]

GROUPS= config["GROUPS"]

CONTRASTS= config["CONTRASTS"]

# Genome, gtf and annotation file for deg 
 
GENOMEFILE= config["GENOMEFILE"]

GTFFILE= config["GTFFILE"]

ANNOTATE= config["ANNOTATE"]

# Running options (yes/no)

DEG= config["DEG"] 

TRIM= config["TRIM"]

# resume wf after qc failure  
RESUME= config["RESUME"]


# Trimmomatic options

FASTAWITHADAPTERSETC= config["FASTAWITHADAPTERSETC"]

SEEDMISMATCHES= config["SEEDMISMATCHES"]

PALINDROMECLIPTHRESHOLD= config["PALINDROMECLIPTHRESHOLD"]

SIMPLECLIPTHRESHOLD= config["SIMPLECLIPTHRESHOLD"]

WINDOWSIZE= config["WINDOWSIZE"]

WINDOWQUALITY= config["WINDOWQUALITY"]

LEADINGQUALITY= config["LEADINGQUALITY"]

TRAILINGQUALITY= config["TRAILINGQUALITY"]

CROPLENGTH= config["CROPLENGTH"]

HEADCROPLENGTH= config["HEADCROPLENGTH"]

MINLEN= config["MINLEN"]

TARGETLENGTH= config["TARGETLENGTH"]

STRICTNESS= config["STRICTNESS"]

# STAR options

STARDIR= config["STARDIR"]

STARVER= config["STARVER"]

SJDBOVERHANG= config["SJDBOVERHANG"]

OUTSAMUNMAPPED= config["OUTSAMUNMAPPED"]

ADAPTER1= config["ADAPTER1"]

ADAPTER2= config["ADAPTER2"]

# cufflinks options (for now defaults)
FILTERINTRONMOTIFS= config["FILTERINTRONMOTIFS"]
SAMSTRANDFIELD = config["SAMSTRANDFIELD"]
# to use as Encode options (for now we have defaults except first) 
FILTERTYPE =  config["FILTERTYPE"]
FILTERMULTIMAPNMAX = config["FILTERMULTIMAPNMAX"]
ALIGNSJOVERHANGMIN = config["ALIGNSJOVERHANGMIN"]
ALIGNSJDBOVERHANGMIN=  config["ALIGNSJDBOVERHANGMIN"]
FILTERMISMATCHNMAX=  config["FILTERMISMATCHNMAX"]
FILTERMISMATCHNOVERLMAX = config["FILTERMISMATCHNOVERLMAX"]
ALIGNINTRONMIN =  config["ALIGNINTRONMIN"]
ALIGNINTRONMAX = config["ALIGNINTRONMAX"]
ALIGNMATESGAPMAX =  config["ALIGNMATESGAPMAX"]
# wig options (for now default)
WIGTYPE = config["WIGTYPE"]
WIGSTRAND = config["WIGSTRAND"]


# Tools options

FASTQCVER= config["FASTQCVER"]

CHKQCVER= config["CHKQCVER"]

TRIMMOMATICVER= config["TRIMMOMATICVER"]

PICARDVER= config["PICARDVER"]

DIR= config["MYFOLDER"]

SUBREADVER= config["SUBREADVER"]

# subread option
STRANDED= config["STRANDED"]

BWAVER=config["BWAVER"]

RNASEQCVER= config["RNASEQCVER"]

RRNALIST= config["RRNALIST"]
# filtering low counts

MINCOUNT= config["MINCOUNT"]
MINSAMPLES= config["MINSAMPLES"]

#  Rules -------------------------------------------------------------------- 

if DEG == "yes" and TRIM == "yes":
  rule all:
     params: batch='-l nodes=1:gpfs'
     input: "postTrimQC","files_to_rnaseqc.txt","STAR_QC","sampletable.txt","deseq2_pca.png","edgeR_prcomp.png","RawCountFile_filtered.txt","Limma_MDS.png"

elif TRIM == "yes" and DEG == "no":
  rule all:
     params: batch='-l nodes=1:gpfs'
     input: "postTrimQC","files_to_rnaseqc.txt","STAR_QC","RawCountFile_filtered.txt"
elif TRIM == "no" and DEG == "yes":
  rule all:
     params: batch='-l nodes=1:gpfs'
     input: "files_to_rnaseqc.txt","STAR_QC","RawCountFile_filtered.txt","sampletable.txt","deseq2_pca.png","edgeR_prcomp.png","Limma_MDS.png"
else:
  rule all:
     params: batch='-l nodes=1:gpfs'
     input: "files_to_rnaseqc.txt","STAR_QC","RawCountFile_filtered.txt"
     ## input: "files_to_rnaseqc.txt","STAR_QC"


if TRIM == "yes":
   rule trimmomatic_pe:
      input: file1= INPUTDIR+"{name}_R1_all.fastq.gz",file2=INPUTDIR+"{name}_R2_all.fastq.gz" 
      output: out11="trim/{name}_R1_all_trim_paired.fastq.gz",out12="trim/{name}_R1_all_trim_unpaired.fastq.gz",out21="trim/{name}_R2_all_trim_paired.fastq.gz",out22="trim/{name}_R2_all_trim_unpaired.fastq.gz"
  ## log: "trim.log" avoid for time saving
      params: batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads:32
      shell:"module load {TRIMMOMATICVER}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE -threads {threads} {input.file1} {input.file2}        {output.out11} {output.out12} {output.out21} {output.out22} ILLUMINACLIP:{FASTAWITHADAPTERSETC}:{SEEDMISMATCHES}:{PALINDROMECLIPTHRESHOLD}:{SIMPLECLIPTHRESHOLD}  LEADING:{LEADINGQUALITY} TRAILING:{TRAILINGQUALITY} SLIDINGWINDOW:{WINDOWSIZE}:{WINDOWQUALITY} MAXINFO:{TARGETLENGTH}:{STRICTNESS} MINLEN:{MINLEN} HEADCROP:{HEADCROPLENGTH} "
## missing CROP:{CROPLENGHT}

   rule fastqc:  
      input: expand("trim/{name}_R1_all_trim_paired.fastq.gz", name=SAMPLES), expand("trim/{name}_R2_all_trim_paired.fastq.gz", name=SAMPLES)  
      output: "postTrimQC"
      priority: 2
      params: batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads: 32
      shell: "mkdir {output}; module load {FASTQCVER} ; fastqc {input} -t {threads} -o {output}"

   rule check:
      input: "postTrimQC"
      output: out1="FastqcSummary.xlsx", out2="fastqc_status.txt"
      params: batch='-l nodes=1:gpfs'
      shell: """
             {CHKQCVER} --dir {input} --outfile {output.out1};
             if [ "$?" == 1 ]; then 
                 echo "FAILURE" > {output.out2}
             else 
                 echo "OK" > {output.out2} 
             fi
             """
   rule decide:
      input: file1="fastqc_status.txt"
      output: out1="fastqc_status_checked.txt"
      params: batch='-l nodes=1:gpfs'
      shell: """
             status=$(<{input.file1})
             if [ $status == "OK" -o {RESUME} == "yes" ]; then
                  if [ {RESUME} == "yes" -a $status == "FAILURE"  ]; then
                      echo "workflow resumed after QC check" > {output.out1}
                  else
                      echo "status OK" > {output.out1}
                  fi 
             else
                  exit 1
             fi
             """

   rule star1p:
      input: file1= "trim/{name}_R1_all_trim_paired.fastq.gz",file2="trim/{name}_R2_all_trim_paired.fastq.gz",file3="fastqc_status_checked.txt"
      output: out1= "{name}.SJ.out.tab",out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: prefix="{name}",batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads: 32
      shell: """
             module load {STARVER}
             STAR --genomeDir {STARDIR} --genomeLoad LoadAndKeep --outFilterIntronMotifs {FILTERINTRONMOTIFS} --outSAMstrandField {SAMSTRANDFIELD}  --outFilterType {FILTERTYPE} --outFilterMultimapNmax {FILTERMULTIMAPNMAX} --alignSJoverhangMin {ALIGNSJOVERHANGMIN} --alignSJDBoverhangMin {ALIGNSJDBOVERHANGMIN}  --outFilterMismatchNmax {FILTERMISMATCHNMAX} --outFilterMismatchNoverLmax {FILTERMISMATCHNOVERLMAX}  --alignIntronMin {ALIGNINTRONMIN} --alignIntronMax {ALIGNINTRONMAX} --alignMatesGapMax {ALIGNMATESGAPMAX} --clip3pAdapterSeq {ADAPTER1} {ADAPTER2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}.
             awk 'BEGIN {{OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";}} {{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' {output.out1}  > {output.out2}
             """  

   rule star2p:
      input: file1= "trim/{name}_R1_all_trim_paired.fastq.gz",file2="trim/{name}_R2_all_trim_paired.fastq.gz",dir="STARINDEX"
      output: "{name}.p2.Aligned.out.sam"
      params: prefix="{name}.p2",batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads:32
      shell:"module load {STARVER}; STAR --genomeDir {input.dir} --genomeLoad LoadAndKeep -outFilterIntronMotifs {FILTERINTRONMOTIFS} --outSAMstrandField {SAMSTRANDFIELD}  --outFilterType {FILTERTYPE} --outFilterMultimapNmax {FILTERMULTIMAPNMAX} --alignSJoverhangMin {ALIGNSJOVERHANGMIN} --alignSJDBoverhangMin {ALIGNSJDBOVERHANGMIN}  --outFilterMismatchNmax {FILTERMISMATCHNMAX} --outFilterMismatchNoverLmax {FILTERMISMATCHNOVERLMAX}  --alignIntronMin {ALIGNINTRONMIN} --alignIntronMax {ALIGNINTRONMAX} --alignMatesGapMax {ALIGNMATESGAPMAX}  --clip3pAdapterSeq {ADAPTER1} {ADAPTER2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMunmapped {OUTSAMUNMAPPED} --outWigType {WIGTYPE} --outWigStrand {WIGSTRAND}"


else:
   rule star1p:      
      input: file1= INPUTDIR+"{name}_R1_all.fastq.gz",file2=INPUTDIR+"{name}_R2_all.fastq.gz" 
      output: out1= "{name}.SJ.out.tab",out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: prefix="{name}",batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads: 32
      shell:"module load {STARVER} ; STAR --genomeDir {STARDIR} --genomeLoad LoadAndKeep -outFilterIntronMotifs {FILTERINTRONMOTIFS} --outSAMstrandField {SAMSTRANDFIELD}  --outFilterType {FILTERTYPE} --outFilterMultimapNmax {FILTERMULTIMAPNMAX} --alignSJoverhangMin {ALIGNSJOVERHANGMIN} --alignSJDBoverhangMin {ALIGNSJDBOVERHANGMIN}  --outFilterMismatchNmax {FILTERMISMATCHNMAX} --outFilterMismatchNoverLmax {FILTERMISMATCHNOVERLMAX}  --alignIntronMin {ALIGNINTRONMIN} --alignIntronMax {ALIGNINTRONMAX} --alignMatesGapMax {ALIGNMATESGAPMAX}  --clip3pAdapterSeq {ADAPTER1} {ADAPTER2}  --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}.; awk 'BEGIN {{OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";}} {{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' {output.out1}  > {output.out2}"
   
   rule star2p:
      input: file1= INPUTDIR+"{name}_R1_all.fastq.gz",file2=INPUTDIR+"{name}_R2_all.fastq.gz",dir="STARINDEX"
      output: "{name}.p2.Aligned.out.sam"
      params: prefix="{name}.p2",batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
      threads:32
      shell:"module load {STARVER}; STAR --genomeDir {input.dir} --genomeLoad LoadAndKeep -outFilterIntronMotifs {FILTERINTRONMOTIFS} --outSAMstrandField {SAMSTRANDFIELD}  --outFilterType {FILTERTYPE} --outFilterMultimapNmax {FILTERMULTIMAPNMAX} --alignSJoverhangMin {ALIGNSJOVERHANGMIN} --alignSJDBoverhangMin {ALIGNSJDBOVERHANGMIN}  --outFilterMismatchNmax {FILTERMISMATCHNMAX} --outFilterMismatchNoverLmax {FILTERMISMATCHNOVERLMAX}  --alignIntronMin {ALIGNINTRONMIN} --alignIntronMax {ALIGNINTRONMAX} --alignMatesGapMax {ALIGNMATESGAPMAX}  --clip3pAdapterSeq {ADAPTER1} {ADAPTER2}  --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat  --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMunmapped {OUTSAMUNMAPPED} --outWigType {WIGTYPE} --outWigStrand {WIGSTRAND}"


rule newindex:
  input: expand("{name}.SJ.out.tab.Pass1.sjdb", name=SAMPLES)
  output: dir="STARINDEX",file="all.SJ.out.tab.Pass1.sjdb"
  params: batch='-l nodes=1:gpfs, mem=250gb, ncpus=32'
  threads: 32
  shell:"mkdir {output.dir}; cat {input} > {output.file}; module load {STARVER}; STAR --runMode genomeGenerate --genomeDir {output.dir} --genomeFastaFiles {GENOMEFILE} --sjdbGTFfile {GTFFILE} --sjdbFileChrStartEnd {output.file} --sjdbOverhang {SJDBOVERHANG} --runThreadN {threads}" 


rule picard:
  input: file1= "{name}.p2.Aligned.out.sam"
  output: outstar1="{name}.star_rg_added.sorted.bam", outstar2="{name}.star_rg_added.sorted.dmark.bam",outstar3="{name}.star.duplic" 
  params: batch='-l nodes=1:gpfs:g24:c16'
  shell: "module load {PICARDVER}; java -Xmx10g  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar INPUT={input.file1} OUTPUT={output.outstar1} SORT_ORDER=coordinate RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; java -Xmx10g -jar $PICARDJARPATH/MarkDuplicates.jar INPUT={output.outstar1} OUTPUT={output.outstar2} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3}"


rule prernaseqc:
  input: expand("{name}.star_rg_added.sorted.dmark.bam", name=SAMPLES)
  output: out1="files_to_rnaseqc.txt"
  priority: 2
  params: batch='-l nodes=1:gpfs'
  run:
        with open(output.out1, "w") as out:
            out.write("Sample ID\tBam file\tNotes\n")
            for f in input:
                out.write("%s\t"  % f)
                out.write("%s\t"  % f)
                out.write("%s\n"  % f)
            out.close()

rule rnaseqc:
  input: "files_to_rnaseqc.txt"
  output: "STAR_QC"
  priority: 2
  params: batch='-l nodes=1:gpfs:g24:c16'
#  shell: "module load {BWAVER}; java -jar {RNASEQCVER} -n 1000 -s {input} -t {GTFFILE} -r {GENOMEFILE}  -o {output}"
  shell: """
         module load {BWAVER}
         var="{RRNALIST}"
         if [  $var == "-" ]; then
                java -jar {RNASEQCVER} -n 1000 -s {input} -t {GTFFILE} -r {GENOMEFILE}  -o {output}
         else
                java -jar {RNASEQCVER} -n 1000 -s {input} -t {GTFFILE} -r {GENOMEFILE} -rRNA {RRNALIST}  -o {output}
         fi
         """

rule subread:
   input:  "{name}.star_rg_added.sorted.dmark.bam"
   output: out="{name}.star.count.info.txt", res="{name}.star.count.txt"
   params: batch='-l nodes=1:gpfs:g24:c16'
   shell: "module load {SUBREADVER}; featureCounts -T 16 -s {STRANDED} -p -t exon -g gene_name -a {GTFFILE} -o {output.out}  {input}; sed '1d' {output.out} | cut -f1,7 > {output.res}"


rule rawfileplot: 
   input: files=expand("{name}.star.count.txt", name=SAMPLES)
   output: "RawCountFile_filtered.txt"
   params: batch='-l nodes=1:gpfs:g24:c16'
   run:
        R("""
        library('reshape') 
        library('ggplot2')
        library('edgeR')
        setwd("{DIR}")
        myfiles=as.character(unlist(strsplit("{input.files}", split=" ")))
        res=read.delim(myfiles[1],header=T)
        colnames(res)[1]="gene"
        colnames(res)[2]=as.character(myfiles[1]) 
        # remove the last 5 statistics lines ... 
        # nr=dim(res)[1]
        # res=res[-c((nr-4):nr),]
        #
        for(i in seq(2, length(myfiles), by = 1))
        {{
        temp=read.delim(myfiles[i],header=T)
        colnames(temp)[1]="gene"
        colnames(temp)[2]=as.character(myfiles[i]) 
        res=merge(res,temp)
        }}
        write.table(as.data.frame(res),file="RawCountFile.txt",sep="\t",row.names=F) 
        #
        mydata=read.delim("RawCountFile.txt",row.names=1)
        val1=as.numeric("{MINCOUNT}")
        val2=as.numeric("{MINSAMPLES}")
        cat(val1," ", val2, "checking..\n",file="check.txt")
        filter <- apply(mydata, 1, function(x) length(x[x>val1])>=val2)
        res=mydata[filter,]
        write.table(as.data.frame(res),file="RawCountFile_filtered.txt",sep="\t",col.names=NA)
        png("HistBeforenormFilter.png")
        df.m <- melt(as.data.frame(res))
        print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
        dev.off() 
        y = DGEList(counts=res)
        ## Normalization TMM ------------------------------------------------------------
        ## method = =c("TMM","RLE","upperquartile","none")
        y <- calcNormFactors(y,method="TMM")
        ndata= cpm(y,log=FALSE,normalized.lib.sizes=TRUE)
        ## save it 
        write.table(ndata,file="CPM_TMM_counts.txt",sep="\t",col.names=NA)
        """)

rule samplecondition:
   input: files=expand("{name}.star.count.txt", name=SAMPLES)
   output: out1= "sampletable.txt"
   params: batch='-l nodes=1:gpfs'
   run:
        with open(output.out1, "w") as out:
            out.write("sampleName\tfileName\tcondition\n")
            i=0
            for f in input.files:
                out.write("%s\t"  % f)
                out.write("%s\t"  % f)
                out.write("%s\n" % GROUPS[i])
                i=i+1
            out.close()

rule deseq2:
  input: file1="sampletable.txt", file2="RawCountFile_filtered.txt"
  ## input: "sampletable.txt"
  output: "deseq2_pca.png"
  params: batch='-V -l nodes=1:gpfs:g24:c16'
  run:
        R("""
        library('DESeq2')
        library('RColorBrewer') 
        library('gplots')
        library('reshape') 
        library('ggplot2')

        #
        setwd("{DIR}")
        sampleinfo=read.delim("{input.file1}")
        sampleFiles=as.character(sampleinfo[,2])
        x = read.delim("{input.file2}",row.names=1)
        ## read annotation file
        ## ann=read.delim("{ANNOTATE}")
        #
        ddsHTSeq<-DESeqDataSetFromMatrix(countData=x,colData=sampleinfo, design=~condition)
        dds<-DESeq(ddsHTSeq)
        ndata=as.data.frame(counts(dds,normalized=TRUE))
        colnames(ndata)=colnames(x)
        write.table(ndata,file="Deseq2_normalized_counts.txt",sep="\t",col.names=NA)
        png("HistDesq2normFilter.png")
        df.m <- melt(as.data.frame(ndata))
        print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
        dev.off() 
        #
        contras=unlist(strsplit("{CONTRASTS}", split=" "))        
        cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontra.txt")
        for(i in seq(1, length(contras), by = 2))
        {{
        res<-results(dds,contrast=c("condition",as.character(contras[i]),as.character(contras[i+1])))
        res<-res[order(res$padj),]
        res1=as.data.frame(res)
        write.table(res1,file=paste("DEG_",contras[i],"_vs_",contras[i+1],".txt",sep=""),sep="\t",col.names=NA) 
        x=res1$log2FoldChange[which(!is.na(res1$log2FoldChange))] 
        plotMA(res,ylim=range(x),main=paste("MAplot_",contras[i],"_vs_",contras[i+1],sep=""))
        dev.copy(png,paste("MAplot_",contras[i],"_vs_",contras[i+1],".png",sep=""))
        dev.off()
        }}
        ## transformation
        rld <- rlogTransformation(dds, blind=TRUE)
        rldm=assay(rld)
        colnames(rldm)=colnames(x)
        write.table(rldm,file="Deseq2_normalized_rld.txt",sep="\t",col.names=NA)
        ## clustering
        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        distsRL <- dist(t(assay(rld)))
        mat <- as.matrix(distsRL)
        rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition,sampleFiles , sep=" : "))
        #if you just want the conditions use this line : rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
        heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16))
        dev.copy(png,"deseq2_heatmaps_samplebysample.png")
        dev.off()
        ## plotMA(dds,ylim=c(-2,2),main="DESeq2 MAplot")
        ## dev.copy(png,"deseq2_MAplot.png")
        ## dev.off()
        ## pca
        print(plotPCA(rld, intgroup=c("condition")))
        dev.copy(png,"deseq2_pca.png")
        dev.off()
        """)

rule edgeR:
  input: file1="sampletable.txt", file2="RawCountFile_filtered.txt"
  output: "edgeR_prcomp.png"
  params: batch='-V -l nodes=1:gpfs:g24:c16'
  run:
        R("""
        library('edgeR')
        library('statmod')
        library('RColorBrewer') 
        library('gplots')
        library('reshape') 
        library('ggplot2')

        #
        setwd("{DIR}")
        # read files
        sampleinfo=read.delim("{input.file1}")
        x = read.delim("{input.file2}",row.names=1)
        # sampleFiles=as.character(sampleinfo[,2])
        ## read annotation file
        ## ann=read.delim("{ANNOTATE}")
        # DGElist object --------------------------------------------------------------
        condition = as.factor(sampleinfo$condition)
        y = DGEList(counts=x,group=condition)
        ## Normalization TMM ------------------------------------------------------------
        ## method = =c("TMM","RLE","upperquartile","none")
        y <- calcNormFactors(y,method="TMM")
        # y$samples
        png("libdistrib.png")
        barplot(y$samples$lib.size*1e-6,main="Library size distribution", names= strsplit(colnames(y$counts),".star.count.txt"), ylab="Library size (millions)",las=2,cex.names=0.8)
        dev.off()
        ## MDS plots ----------------------------------------------------------------------
        # both pairewise (leading)
        png("MDS_bcv.png")
        plotMDS(y, method="bcv", , main="MDS plot bcv")
        dev.off()
        png("MDS_logFC.png")
        plotMDS(y, method="logFC" , main="MDS plot logFC") ## plotMDS(y) default
        dev.off()
        # plotMDS(y, method="logFC",gene.selection="common", main="MDS plot common")
        ## estimating common and tagwise dispersions -----------------------------------------
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y) #default trend: moveingave
        ## plotting
        png("BCVplot.png")
        plotBCV(y,main="BCV plot")
        dev.off()
        ## differentially expressed genes ---------------------------------------------------
        contras=unlist(strsplit("{CONTRASTS}", split=" "))        
        cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontra-edgeR.txt")
        for(i in seq(1, length(contras), by = 2))
        {{
        deg<-exactTest(y,c(as.character(contras[i+1]),as.character(contras[i])))
        # 
        n=dim(y$counts)[1]
        tt=topTags(deg, n=n)
        res1 = as.data.frame(tt)
        #
        ## res1=cbind(Ensembl.Gene.ID=substr(rownames(res1),1,18),id.ver=rownames(res1),res1)
        ## final= merge(res1,ann,all.x=TRUE)
        ##final=final[order(final$FDR),]
        final=res1[order(res1$FDR),]
        write.table(final,file=paste("DEG_EdgeR_",contras[i],"_vs_",contras[i+1],".txt",sep=""),sep="\t",col.names=NA)
        #  like MAplot
        deg1sel <- decideTestsDGE(deg, p=0.05, adjust="BH")
        detags <- rownames(y)[as.logical(deg1sel)]
        png(paste("Smearplot_",contras[i],"_vs_",contras[i+1],".png",sep=""))
        plotSmear(deg, de.tags=detags,main= paste("Smearplot FDR<0.05 ",contras[i],"_vs_",contras[i+1],sep=""))
        abline(h = c(-2, 2), col = "blue")
        dev.off()
        # 
        }}
        ## transformation
        ylog2=cpm(y,log=TRUE,normalized.lib.sizes=TRUE,prior.count=2) # prior count like avelogcpm
        ndata= cpm(y,log=FALSE,normalized.lib.sizes=TRUE)*1e6
        ## save it
        write.table(ylog2,file="edgeR_normalized_counts_log.txt",sep="\t",col.names=NA) 
        write.table(ndata,file="edgeR_normalized_counts.txt",sep="\t",col.names=NA)
        png("HistEdgeRnormFilter.png")
        df.m <- melt(as.data.frame(ndata))
        print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
        dev.off()         
        ## clustering / heatmap
        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        distylog2=dist(t(ylog2))
        mat = as.matrix(distylog2)
        # rownames(mat) <- colnames(mat)
        heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16))
        dev.copy(png,"edgeR_heatmaps_samplebysample.png")
        dev.off()
        #pca
        pr2=prcomp(t(ylog2))
        png("edgeR_prcomp.png")
        # biplot(pr2)
        plot(pr2$x[,1],pr2$x[,2],col="red", main="PCA plot using prcomp and Logcpm data")
        text(pr2$x[,1],pr2$x[,2], labels=colnames(ylog2), cex=0.7, pos=4)
        dev.off()
        """)

rule limmavoom:
  input: file1="sampletable.txt", file2="RawCountFile_filtered.txt"
  output: "Limma_MDS.png"
  params: batch='-V -l nodes=1:gpfs:g24:c16'
  run:
        R("""
        library('edgeR')
        library('statmod')
        library('RColorBrewer') 
        library('gplots')
        library('reshape') 
        library('ggplot2')
        library('limma')
        library('geneplotter')
        setwd("{DIR}")
        # read files
        sampleinfo=read.delim("{input.file1}")
        x = read.delim("{input.file2}",row.names=1)
        # sampleFiles=as.character(sampleinfo[,2])
        Group <- factor(sampleinfo$condition)
        design=model.matrix(~0+Group)
        contras=unlist(strsplit("{CONTRASTS}", split=" "))        
        cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontraLimma.txt")
        cons=c()
        for(i in seq(1, length(contras), by = 2))
        {{
        cons=c(cons,paste(contras[i],"-",contras[i+1],sep=""))
        }}
        v1 <- voom(as.matrix(x),design,plot=TRUE,normalize="quantile")
        sf = v1$E/log2((x/colSums(x))*1000000)
        write.table(sf,file="LimmaVoom_scaling_factors.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
        png("HistLimmavoomNormFilter.png")
        df.n <- melt(as.data.frame(v1$E))
        print(ggplot(df.n) + geom_density(aes(x = value,colour = variable)) + labs(x = NULL) + theme(legend.position='top'))
        dev.off()
        ## MDS
        MDS <- plotMDS(v1,xlim=c(-5,5),ylim=c(-5,5),cex=1,pch=20)
        shortname=paste(substr(colnames(v1$E),1,22))
        text(MDS, labels=shortname, cex=0.5, pos=1)
        print(MDS)
        dev.copy(png, paste("Limma_MDS.png"))
        dev.off()
        ## DEG
        nb=length(contras)/2
        colnames(design) <- levels(Group)
        fit <- lmFit(v1,design)
        contrast.matrix <- makeContrasts(contrasts=cons,levels=design)
        fitb <- contrasts.fit(fit, contrast.matrix)
        ebayes.fit=eBayes(fitb)
        ## 
        for (i in 1:nb)
        {{
        all.genes.con = topTable(ebayes.fit, coef = i, number=nrow(ebayes.fit))
        ## generate Volcano plot
        jpeg(paste("Limma_",cons[i],"_volcano.jpeg",sep=""),quality=100) 
        plot(all.genes.con$logFC,-log10(all.genes.con$adj.P.Val),cex=0.1,xlab="Log Fold-Change",ylab="-log10 Adj P-Value",main=paste('Volcano Plot for ',cons[i],sep=""))
        t=which(all.genes.con$adj.P.Val<0.05 & abs(all.genes.con$logFC)>=1 )
        points(all.genes.con$logFC[t],-log10(all.genes.con$adj.P.Val[t]),col="red",pch=20,cex=0.5)
        dev.off()
        #MAplot <- plot(ebayes.fit,coef=i)
        #print(MAplot)
        #dev.copy(png, paste(cons[i],"_MAplot_Limma_old.png",sep=""))
        #dev.off()
        dataf=data.frame("m"=all.genes.con$AveExpr,"fc"=all.genes.con$logFC,"sig"=all.genes.con$adj.P.Val<0.05)
        png(paste(cons[i],"_MAplot_Limma_v2.png",sep=""))
        plotMA(dataf,log="",main=cons[i],ylim=range(all.genes.con$logFC))
        dev.off()
        write.table(all.genes.con,file=paste("Limma_deg_",cons[i],"_all_genes.txt",sep=""),sep="\t",col.names=NA)
        }}
        #
        """)


