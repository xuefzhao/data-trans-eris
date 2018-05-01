#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-p", "--path"), type="character", default=NULL,
              help="path of input bed files", metavar="character"),
        make_option(c("-o","--output"), type="character", default=NULL,
              help="name of output stat file", metavar="character" )
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out=data.frame('fam'=0,'ro'=0,'denovo'=0,'fa_pb'=0,'mo_pb'=0,'fa_mo_pb'=0)
out2=data.frame('SVTYPE'=0,'fam'=0,'ro'=0,'denovo'=0,'fa_pb'=0,'mo_pb'=0,'fa_mo_pb'=0)
rec=0
dat_svtype_modify<-function(dat){
	dat[,ncol(dat)+1]=0
	for(i in unique(dat[,4])){
		i_new=strsplit(i,'_')[[1]][1]
		if(i=='ALU' | i=='LINE1' | i=='SVA'){	i_new='INS'}
		dat[dat[,4]==i,][,ncol(dat)]=i_new
	}
	dat[,4]=dat[,ncol(dat)]
	return(dat[,c(1:(ncol(dat)-1))])
	}
for(k1 in list.files(opt$path)){
        if(length(strsplit(k1,'[.]')[[1]])>1 & strsplit(k1,'[.]')[[1]][2]=='pb'){
                dat=read.table(paste(opt$path,k1,sep='/'))
                dat=dat_svtype_modify(dat)
                stat=data.frame(table(dat[,ncol(dat)]))
                rec=rec+1
                out[rec,1]=strsplit(k1,'[.]')[[1]][1]
                out[rec,2]=strsplit(k1,'[.]')[[1]][3]
                out[rec,3]=stat[stat[,1]=='denovo',2]
                out[rec,4]=stat[stat[,1]=='fa_pb',2]
                out[rec,5]=stat[stat[,1]=='mo_pb',2]
                out[rec,6]=stat[stat[,1]=='fa_mo_pb',2]
                stat2=cbind(strsplit(k1,'[.]')[[1]][1],strsplit(k1,'[.]')[[1]][3],table(dat[,c(4,ncol(dat))]))
                stat2=cbind(rownames(stat2),stat2)
                colnames(stat2)[c(1:3)]=c('SVTYPE','fam','ro')
                out2=rbind(out2, stat2)
        }
}

out=out[order(out[,2]),]
out=out[order(out[,1]),]
write.table(out,opt$output, quote=F, sep='\t', col.names=T, row.names=F)
write.table(out2[-1,],paste(opt$output,'.by_svtype',sep=''), quote=F, sep='\t', col.names=T, row.names=F)

