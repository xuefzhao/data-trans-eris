#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-f", "--father"), type="character", default=NULL,
              help="name of father bed", metavar="character"),
        make_option(c("-m", "--mother"), type="character", default=NULL,
              help="name of mother bed", metavar="character"),
        make_option(c("-p", "--proband"), type="character", default=NULL,
              help="name of proband bed", metavar="character"),
        make_option(c("-o","--output"), type="character", default=NULL,
              help="name of output file", metavar="character" )
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fa=read.table(opt$father)
mo=read.table(opt$mother)
pb=read.table(opt$proband)
fa[,ncol(fa)+1]=paste(fa[,1],fa[,2],fa[,3],fa[,4],sep='_')
mo[,ncol(mo)+1]=paste(mo[,1],mo[,2],mo[,3],mo[,4],sep='_')
pb[,ncol(pb)+1]='denovo'
pb[,ncol(pb)+1]=paste(pb[,1],pb[,2],pb[,3],pb[,4],sep='_')

pb[pb[,ncol(pb)]%in%fa[,ncol(fa)],][,ncol(pb)-1]='fa_pb'
pb[pb[,ncol(pb)]%in%mo[,ncol(mo)],][,ncol(pb)-1]='mo_pb'
pb[pb[,ncol(pb)]%in%fa[,ncol(fa)] & pb[,ncol(pb)]%in%mo[,ncol(mo)],][,ncol(pb)-1]='fa_mo_pb'

write.table(pb[,c(1:(ncol(pb)-1))], opt$output, quote=F, sep='\t', col.names=F, row.names=F)


