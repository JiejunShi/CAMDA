#### libraries
library("getopt");
#### arguments
#list of options:long name, short name, flag(0=no argument, 1=required argument, 2=optional argument), data type, description
OptList=c('help',       'h', 0, 'logical',   'useage',
          'ReadCT',     'i', 1, 'character', 'ReadCT file. REQUIRED.',
          'Regions',    'r', 1, 'character', 'Bed file of regions whose Methylation Concurrence Ratio(CAMDA) will be reported. REQUIRED.',
          'UseStrand',  's', 0, 'logical', 'If -s is specified, strand infomation(6th column) in Regions file will be used.',
          'Weight',     'w', 2, 'character', 'Weight applied to each sub-read fragment, either "cg" or "1". "cg" means weight equal to the CpG count of each fragment. "1" means no weight applied. [Default="cg"]',
          'Output',     'o', 2, 'character', 'Output file report CAMDA and MethRatio of each region. [Default="Region_CAMDA.tsv"].'
        );
OptList=matrix(OptList, byrow=TRUE, ncol=5)
opt=getopt(OptList);
################
if(!is.null(opt$help)){
  cat(getopt(OptList, command = paste("Rscript",get_Rscript_filename()), usage=TRUE));
  q(status=1);
};
if(is.null(opt$ReadCT)){
  print("[ERROR] No input ReadCT file detected.");
  q(status=1)
}else{
  ReadCT=as.character(opt$ReadCT)
};
if(is.null(opt$Regions)){
  print("[ERROR] No input Regions bed file detected.");
  q(status=1)
}else{
  Regions=as.character(opt$Regions)
};
if(!is.null(opt$UseStrand)){UseStrand=TRUE}else{UseStrand=FALSE};
if(is.null(opt$Weight)){
  Weight="cg";
}else{
  Weight=as.character(opt$Weight);
  if(Weight!="cg" & Weight!="1"){
    print("[ERROR] Weight should be either cg or 1.");
    q(status=1)
  };
};
if(is.null(opt$Output)){
  Output="Region_CAMDA.tsv";
}else{
  Output=as.character(opt$Output)
};
################
# consecutive C/T sub-read fragment count
CTmat2MR_CAMDA<-function(x,w="cg"){
  methyl_frag=c();unmethyl_frag=c();interrupt_frag=c();
  for(i in 1:length(x)){
    C_frag=attributes(gregexpr("C+",x[i])[[1]])$match.length;
    #if(length(C_frag)==1 & C_frag==(-1)){C_frag=0};
    T_frag=attributes(gregexpr("T+",x[i])[[1]])$match.length;
    #if(length(T_frag)==1 & T_frag==(-1)){T_frag=0};
    I_frag=attributes(gregexpr("I+",x[i])[[1]])$match.length;
    if(sum(C_frag)==(-1) & sum(T_frag)==(-1) & sum(I_frag)==(-1)){
      next
    }else if(sum(C_frag)==(-1) & sum(T_frag)>0){
      unmethyl_frag=c(unmethyl_frag,T_frag)
    }else if(sum(C_frag)>0 & sum(T_frag)==(-1)){
      methyl_frag=c(methyl_frag,C_frag)
    }else if(sum(C_frag)>0 & sum(T_frag)>0){
      methyl_frag=c(methyl_frag,C_frag)
      interrupt_frag=c(interrupt_frag,T_frag)
    }else if(sum(I_frag)>0){
      interrupt_frag=c(interrupt_frag,I_frag)
    };
  };
  MethRatio=sum(methyl_frag)/(sum(methyl_frag)+sum(unmethyl_frag)+sum(interrupt_frag));
  if(w=="cg"){
    CAMDA=sum(interrupt_frag)/(sum(methyl_frag)+sum(unmethyl_frag)+sum(interrupt_frag));
  }else if(w=="1"){
    CAMDA=length(interrupt_frag)/(length(methyl_frag)+length(unmethyl_frag)+length(interrupt_frag));
  };
  return(c(MethRatio,CAMDA));
};
################
print(paste("[",Sys.time(),"] Starting"))

if(UseStrand==FALSE){
  Regions=read.table(Regions,sep="\t",header=F)[,1:3];
  Regions=dplyr::distinct(Regions);
  names(Regions)=c("chr","start","end");
}else{
  Regions=read.table(Regions,sep="\t",header=F)[,c(1:3,6)];
  Regions=dplyr::distinct(Regions);
  names(Regions)=c("chr","start","end","strand");
};
Regions=Regions[order(Regions$chr,Regions$start),];
print(paste("[",Sys.time(),"]",nrow(Regions),"Regions imported"))

if(grepl("*.tsv.gz$",ReadCT)==T){
  ReadCT=read.table(gzfile(ReadCT),sep="\t",header=F);
}else{
  ReadCT=read.table(ReadCT,sep="\t",header=F);
};
names(ReadCT)=c("chr","strand","FirstCT","LastCT","CT_count","CT_pos","CT_seq");
ReadCT=ReadCT[ReadCT$CT_count>0,];
print(paste("[",Sys.time(),"] ReadCT file imported"))

if(UseStrand==FALSE){
  write("chr\tstart\tend\tRead_count\tCT_count\tMethRatio\tCAMDA",file = Output)
}else{
  write("chr\tstart\tend\tstrand\tRead_count\tCT_count\tMethRatio\tCAMDA",file = Output)
};
for(i in 1:nrow(Regions)){
  chr0=as.character(Regions$chr[i]);start0=as.numeric(Regions$start[i]);end0=as.numeric(Regions$end[i]);
  if(UseStrand==TRUE){
    strand0=as.character(Regions$strand[i]);
    i1=which(ReadCT$chr==chr0 & ReadCT$strand==strand0 & ((ReadCT$FirstCT<=start0 & ReadCT$LastCT>=start0) | (ReadCT$FirstCT>=start0 & ReadCT$FirstCT<=end0)));
  }else{
    i1=which(ReadCT$chr==chr0 & ((ReadCT$FirstCT<=start0 & ReadCT$LastCT>=start0) | (ReadCT$FirstCT>=start0 & ReadCT$FirstCT<=end0)));
  };
  read_count=length(i1);
  #if(read_count>1e5){i1=sample(i1,1e5)};
  if(read_count>0){
    ctmat=c();region_CT_count=c();
    for(j in 1:length(i1)){
      CT_pos=as.numeric(unlist(strsplit(as.character(ReadCT$CT_pos[i1[j]]),",")));
      region_CT_count=c(region_CT_count,CT_pos);
      CT_seq=unlist(strsplit(as.character(ReadCT$CT_seq[i1[j]]),","));
      inside=CT_seq[which(CT_pos>=start0 & CT_pos<=end0)];
      outside=CT_seq[which(CT_pos<start0 | CT_pos>end0)];
      if(!("C"%in%inside) & ("C"%in%outside)){
        inside=rep("I",length(inside));
      };
      CT_seq=paste(inside,collapse="");
      if(CT_seq!=""){ctmat=c(ctmat,CT_seq)};
    };
    ratios=CTmat2MR_CAMDA(x=ctmat,w=Weight);
    region_CT_count=unique(region_CT_count);region_CT_count=length(which(region_CT_count>=start0 & region_CT_count<=end0));
  }else{
    ratios=c("NA","NA");region_CT_count=0;
  };
  if(UseStrand==TRUE){
    aRecord=c(chr0,start0,end0,strand0,read_count,region_CT_count,ratios);
  }else{
    aRecord=c(chr0,start0,end0,read_count,region_CT_count,ratios);
  };
  write(paste(aRecord,sep="",collapse = "\t"),file = Output,append = T);
};

print(paste("[",Sys.time(),"] Finished"));
