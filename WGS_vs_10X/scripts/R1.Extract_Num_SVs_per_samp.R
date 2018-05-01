#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-p", "--bedpath"), type="character", default=NULL,
              help="path of input bed file", metavar="character"),
        make_option(c("-a", "--appdix"), type="character", default=NULL,
              help="appdix of bed files", metavar="character"),
        make_option(c("-o","--output"), type="character", default=NULL,
              help="name of output", metavar="character" )

 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

bedpath=opt$bedpath
appdix=opt$appdix

file_names=c()
samp_names=c()
for(i in list.files(bedpath)){
	if(grepl(appdix, i)){
		sample_name=(gsub(appdix,'',i))
		samp_names=c(samp_names,sample_name)
		file_names=c(file_names, paste(bedpath,i,sep='/'))		
	}
}


svtype=c()
for(i in file_names){
	dat=read.table(i)
	svtype=unique(dat[,4])
}


out_stat=data.frame(matrix(ncol=length(svtype), nrow=1))
colnames(out_stat)=svtype
for(i in 1:length(file_names)){
	file=file_names[i]
	samp=samp_names[i]
	dat=read.table(file)	
	dat=dat[dat[,3]-dat[,2]>49,]
	stat=data.frame(table(dat[,4]))
	for(j in stat[,1]){
			out_stat[i,colnames(out_stat)==j]=stat[stat[,1]==j,2]
	}
	rownames(out_stat)[i]=samp
}

write.table(out_stat, file=opt$output, quote=F, sep='\t', col.names=T, row.names=T)

