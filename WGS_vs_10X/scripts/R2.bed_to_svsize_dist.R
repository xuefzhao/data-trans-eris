#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-b", "--bedname"), type="character", default=NULL,
              help="name of input bed file", metavar="character"),
        make_option(c("-o", "--outname"), type="character", default=NULL,
              help="name of output file", metavar="character")

 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dat_to_sv_dist<-function(tmp){
        out=data.frame('size'=0,'num'=0)
        seg=0.1
        seg_list=seq(0,6,by=seg)
        tmp[,ncol(tmp)+1]=abs(tmp[,3]-tmp[,2])
        tmp=tmp[tmp[,ncol(tmp)]>0,]
        tmp[,ncol(tmp)+1]=log10(tmp[,ncol(tmp)])
        rec=0
        for(i in seg_list[2:length(seg_list)]){
                tmp2=tmp[tmp[,ncol(tmp)]< i & !tmp[,ncol(tmp)]< i-seg ,]
                rec=rec+1
                out[rec,1]=i - seg/2
                out[rec,2]=nrow(tmp2)
        }
        return(out)
        }

dat_to_sv_dist_overall<-function(dat){
        all_stats=data.frame('size'=0,'num'=0,'svtype'=0)
        dat[,ncol(dat)+1]=0
        for(i in unique(dat[,4])){
            dat[dat[,4]==i,][,ncol(dat)]=strsplit(i,'_')[[1]][1]
            }
        dat[,4]=dat[,ncol(dat)]
        for(i in unique(dat[,4])){
                tmp=dat[dat[,4]==i,]
                stat1=dat_to_sv_dist(tmp)
                stat1[,3]=i
                colnames(stat1)[3]='svtype'
                all_stats=rbind(all_stats, stat1)
                }
        tmp=dat[dat[,4]%in%c('DUP','INS','ALU','LINE1','SVA','HERVK'),]
        stat1=dat_to_sv_dist(tmp)
        stat1[,3]='DUP_INS'
        colnames(stat1)[3]='svtype'
        all_stats=rbind(all_stats, stat1)
        return(all_stats)
        }

file_in=opt$bedname
fileout=opt$outname
dat=read.table(file_in)
svsize_stat=dat_to_sv_dist_overall(dat)
colnames(svsize_stat)[1]='#size'
write.table(svsize_stat, fileout, quote=F, sep='\t', col.names=T, row.names=F)

