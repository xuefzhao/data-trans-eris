order_bed<-function(gatk_chs){
	gatk_chs=gatk_chs[order(gatk_chs[,3]),]
	gatk_chs=gatk_chs[order(gatk_chs[,1]),]
	gatk_chs=gatk_chs[order(gatk_chs[,1]),]
	return(gatk_chs)
	}
add_callset_supp<-function(gatk_chs, callset, sample,dat_path){
	supp1=read.table(paste(dat_path, '/callset_compare/',callset,'.',sample,'.RO00.bed',sep=''))
	supp2=read.table(paste(dat_path, '/callset_compare/',callset,'.',sample,'.RO01.bed',sep=''))
	supp3=read.table(paste(dat_path, '/callset_compare/',callset,'.',sample,'.RO05.bed',sep=''))
	supp1[,ncol(supp1)+1]=paste(supp1[,1],supp1[,2],supp1[,3],supp1[,4],sep='_')
	supp2[,ncol(supp2)+1]=paste(supp2[,1],supp2[,2],supp2[,3],supp2[,4],sep='_')
	supp3[,ncol(supp3)+1]=paste(supp3[,1],supp3[,2],supp3[,3],supp3[,4],sep='_')
	gatk_chs[,ncol(gatk_chs)+1]=-1
	gatk_chs[gatk_chs$SVID%in%supp1[,ncol(supp1)],][,ncol(gatk_chs)]='0-.1'
	gatk_chs[gatk_chs$SVID%in%supp2[,ncol(supp2)],][,ncol(gatk_chs)]='.1-.5'
	gatk_chs[gatk_chs$SVID%in%supp3[,ncol(supp3)],][,ncol(gatk_chs)]='.5-1'
	colnames(gatk_chs)[ncol(gatk_chs)]=paste(strsplit(callset,'_')[[1]][2],'supp',sep='_')
	return(gatk_chs)
	}
add_PB_Vali_to_GATK<-function(gatk_chs,callset, sample, seq){
	PB_vali_raw=read.table(paste(s04_path,'/',callset,'.',sample,'.PB_Vali.',seq,'.bed.vapor',sep=''), header=T)
	PB_vali_raw=unique(order_bed(PB_vali_raw))

	gatk_chs[,ncol(gatk_chs)+1]=paste(gatk_chs[,1],gatk_chs[,2],gatk_chs[,3],sep='_')
	PB_vali_raw[,ncol(PB_vali_raw)+1]=paste(PB_vali_raw[,1],PB_vali_raw[,2],PB_vali_raw[,3],sep='_')
	PB_vali_raw[,4]=as.character(PB_vali_raw[,4])
	PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]%in%gatk_chs[,ncol(gatk_chs)],][,4]=as.character(gatk_chs[gatk_chs[,ncol(gatk_chs)]%in%PB_vali_raw[,ncol(PB_vali_raw)],][,4])
	gatk_chs=gatk_chs[,c(1:(ncol(gatk_chs)-1))]
	PB_vali_raw=PB_vali_raw[,c(1:(ncol(PB_vali_raw)-1))]

	PB_vali_raw[,ncol(PB_vali_raw)+1]=paste(PB_vali_raw[,1],PB_vali_raw[,2],PB_vali_raw[,3],PB_vali_raw[,4],sep='_')
	gatk_chs[,ncol(gatk_chs)+1]='-1'
	gatk_chs[gatk_chs$SVID%in%PB_vali_raw[,ncol(PB_vali_raw)] ,][,ncol(gatk_chs)]=as.character(PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]%in%gatk_chs$SVID,]$VaPoR_GT)
	colnames(gatk_chs)[ncol(gatk_chs)]=paste('PB_vali_',seq,sep='')
	return(gatk_chs)
	}
add_PB_Vali_to_PB<-function(pb_dat,callset, sample, seq){
	PB_vali_raw=read.table(paste(s04_path,'/',callset,'.',sample,'.PB_Vali.',seq,'.bed.vapor',sep=''), header=T)
	PB_vali_raw=unique(order_bed(PB_vali_raw))

	pb_dat[,ncol(pb_dat)+1]=paste(pb_dat[,1],pb_dat[,2],pb_dat[,3],sep='_')
	PB_vali_raw[,ncol(PB_vali_raw)+1]=paste(PB_vali_raw[,1],PB_vali_raw[,2],PB_vali_raw[,3],sep='_')
	qc_tmp=data.frame(table(PB_vali_raw[,ncol(PB_vali_raw)]))[data.frame(table(PB_vali_raw[,ncol(PB_vali_raw)]))[,2]>1,]
	for(i in qc_tmp[,1]){
		for(j in rownames(PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]==i, ])[2:length(rownames(PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]==i, ]))]){
			PB_vali_raw=PB_vali_raw[rownames(PB_vali_raw)!=j,]
		}
	}

	PB_vali_raw[,4]=as.character(PB_vali_raw[,4])
	PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]%in%pb_dat[,ncol(pb_dat)],][,4]=as.character(pb_dat[pb_dat[,ncol(pb_dat)]%in%PB_vali_raw[,ncol(PB_vali_raw)],][,5])
	pb_dat=pb_dat[,c(1:(ncol(pb_dat)-1))]
	PB_vali_raw=PB_vali_raw[,c(1:(ncol(PB_vali_raw)-1))]

	PB_vali_raw[,ncol(PB_vali_raw)+1]=paste(PB_vali_raw[,1],PB_vali_raw[,2],PB_vali_raw[,3],PB_vali_raw[,4],sep='_')
	pb_dat[,ncol(pb_dat)+1]='-1'
	pb_dat[pb_dat$SVID%in%PB_vali_raw[,ncol(PB_vali_raw)] ,][,ncol(pb_dat)]=as.character(PB_vali_raw[PB_vali_raw[,ncol(PB_vali_raw)]%in%pb_dat$SVID,]$VaPoR_GT)
	colnames(pb_dat)[ncol(pb_dat)]=paste('PB_vali_',seq,sep='')
	return(pb_dat)
	}
add_SV_Size<-function(gatk_chs){
	gatk_chs[,ncol(gatk_chs)+1]='over30Kb'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<30000,][,ncol(gatk_chs)]='10Kb-30Kb'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<10000,][,ncol(gatk_chs)]='5Kb-10Kb'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<5000,][,ncol(gatk_chs)]='1Kb-5Kb'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<1000,][,ncol(gatk_chs)]='500bp-1Kb'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<500,][,ncol(gatk_chs)]='100bp-500bp'
	gatk_chs[gatk_chs[,3]-gatk_chs[,2]<100,][,ncol(gatk_chs)]='50bp-100bp'
	if (nrow(gatk_chs[gatk_chs[,3]-gatk_chs[,2]<50,])>0){
		gatk_chs[gatk_chs[,3]-gatk_chs[,2]<50,][,ncol(gatk_chs)]='under50bp'}
	colnames(gatk_chs)[ncol(gatk_chs)]='SVSize_Cate'
	return(gatk_chs)
	}
label_GATK<-function(sample,fam){
	gatk_chs=read.table(paste(dat_path,'/GATK.',sample,'.bed',sep=''))
	colnames(gatk_chs)=c('CHR','POS','END','SVTYPE','SVTYPE2','END2','SVLEN','BP_Accu','Homo_Len')
	gatk_chs=order_bed(gatk_chs)
	gatk_chs=add_SV_Size(gatk_chs)

	inheri=read.table(paste(dat_path,'/GATK.',fam,'.pb.RO05.bed',sep=''))
	inheri=order_bed(inheri)
	gatk_chs[,ncol(gatk_chs)+1]=inheri[,ncol(inheri)]
	colnames(gatk_chs)[ncol(gatk_chs)]='inheri'

	gatk_chs[,ncol(gatk_chs)+1]=paste(gatk_chs[,1],gatk_chs[,2],gatk_chs[,3],gatk_chs[,4],sep='_')
	colnames(gatk_chs)[ncol(gatk_chs)]='SVID'

	gatk_chs=add_callset_supp(gatk_chs, 'GATK_HGSV',sample, dat_path)
	gatk_chs=add_callset_supp(gatk_chs, 'GATK_SVpipe',sample, dat_path)
	gatk_chs=add_callset_supp(gatk_chs, 'GATK_PBNew',sample, dat_path)

	gc=read.table(paste(dat_path,'/genomic_context/GATK.',sample,'.anno.bed',sep=''))
	gc=unique(order_bed(gc))
	gc[,ncol(gc)+1]=paste(gc[,1],gc[,2],gc[,3],gc[,4],sep='_')
	gatk_chs[,ncol(gatk_chs)+1]='-1'
	colnames(gatk_chs)[ncol(gatk_chs)]='GC'
	gatk_chs=gatk_chs[,-9]
	gatk_chs=unique(gatk_chs)
	gatk_chs[gatk_chs$SVID %in% gc[,ncol(gc)],][,ncol(gatk_chs)]=gc[,ncol(gc)-1]

	gatk_chs=add_PB_Vali_to_GATK(gatk_chs,'GATK',sample,'raw')
	gatk_chs=add_PB_Vali_to_GATK(gatk_chs,'GATK',sample,'h0')
	gatk_chs=add_PB_Vali_to_GATK(gatk_chs,'GATK',sample,'h1')

	write.table(gatk_chs,paste(dat_path,'/GATK.',sample,'.labeled.bed',sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}
label_HGSV<-function(sample,fam){
	hgsv_dat=read.table(paste(dat_path,'/HGSV.',sample,'.bed',sep=''))
	colnames(hgsv_dat)=c('CHR','POS','END','SVTYPE','SVLEN','FILTER','CALLER','NumCaller','Other','SVID','SVTYPE2')
	hgsv_dat=order_bed(hgsv_dat)
	hgsv_dat=add_SV_Size(hgsv_dat)

	inheri=read.table(paste(dat_path,'/HGSV.',fam,'.pb.RO05.bed',sep=''))
	inheri=order_bed(inheri)
	hgsv_dat[,ncol(hgsv_dat)+1]=inheri[,ncol(inheri)]
	colnames(hgsv_dat)[ncol(hgsv_dat)]='inheri'

	hgsv_dat$SVID=paste(hgsv_dat[,1],hgsv_dat[,2],hgsv_dat[,3],hgsv_dat[,4],sep='_')

	hgsv_dat=add_callset_supp(hgsv_dat, 'HGSV_GATK',sample, dat_path)
	hgsv_dat=add_callset_supp(hgsv_dat, 'HGSV_SVpipe',sample, dat_path)
	hgsv_dat=add_callset_supp(hgsv_dat, 'HGSV_PBNew',sample, dat_path)

	gc=read.table(paste(dat_path,'/genomic_context/HGSV.',sample,'.anno.bed',sep=''))
	gc=unique(order_bed(gc))
	gc[,ncol(gc)+1]=paste(gc[,1],gc[,2],gc[,3],gc[,4],sep='_')
	hgsv_dat[,ncol(hgsv_dat)+1]='-1'
	colnames(hgsv_dat)[ncol(hgsv_dat)]='GC'
	hgsv_dat=hgsv_dat[,-9]
	hgsv_dat=unique(hgsv_dat)
	hgsv_dat[hgsv_dat$SVID %in% gc[,ncol(gc)],][,ncol(hgsv_dat)]=as.character(gc[,ncol(gc)-1])

	hgsv_dat=add_PB_Vali_to_GATK(hgsv_dat,'HGSV_ILL',sample,'raw')
	hgsv_dat=add_PB_Vali_to_GATK(hgsv_dat,'HGSV_ILL',sample,'h0')
	hgsv_dat=add_PB_Vali_to_GATK(hgsv_dat,'HGSV_ILL',sample,'h1')

	write.table(hgsv_dat,paste(dat_path,'/HGSV.',sample,'.labeled.bed',sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}
label_SVpipe<-function(sample,fam){
	svpipe_dat=read.table(paste(dat_path,'/SVpipe.',sample,'.bed',sep=''))
	colnames(svpipe_dat)=c('CHR','POS','END','SVTYPE','SVLEN','FILTER','CALLER','NumCaller','Other','SVID','SVTYPE2')
	svpipe_dat=order_bed(svpipe_dat)
	svpipe_dat=add_SV_Size(svpipe_dat)

	inheri=read.table(paste(dat_path,'/SVpipe.',fam,'.pb.RO05.bed',sep=''))
	inheri=order_bed(inheri)
	svpipe_dat[,ncol(svpipe_dat)+1]=inheri[,ncol(inheri)]
	colnames(svpipe_dat)[ncol(svpipe_dat)]='inheri'

	svpipe_dat$SVID=paste(svpipe_dat[,1],svpipe_dat[,2],svpipe_dat[,3],svpipe_dat[,4],sep='_')

	svpipe_dat=add_callset_supp(svpipe_dat, 'SVpipe_GATK',sample, dat_path)
	svpipe_dat=add_callset_supp(svpipe_dat, 'SVpipe_HGSV',sample, dat_path)
	svpipe_dat=add_callset_supp(svpipe_dat, 'SVpipe_PBNew',sample, dat_path)

	gc=read.table(paste(dat_path,'/genomic_context/SVpipe.',sample,'.anno.bed',sep=''))
	gc=unique(order_bed(gc))
	gc[,ncol(gc)+1]=paste(gc[,1],gc[,2],gc[,3],gc[,4],sep='_')
	svpipe_dat[,ncol(svpipe_dat)+1]='-1'
	colnames(svpipe_dat)[ncol(svpipe_dat)]='GC'
	svpipe_dat=svpipe_dat[,-9]
	svpipe_dat=unique(svpipe_dat)
	svpipe_dat[svpipe_dat$SVID %in% gc[,ncol(gc)],][,ncol(svpipe_dat)]=as.character(gc[,ncol(gc)-1])

	svpipe_dat=add_PB_Vali_to_GATK(svpipe_dat,'SVpipe',sample,'raw')
	svpipe_dat=add_PB_Vali_to_GATK(svpipe_dat,'SVpipe',sample,'h0')
	svpipe_dat=add_PB_Vali_to_GATK(svpipe_dat,'SVpipe',sample,'h1')

	write.table(svpipe_dat,paste(dat_path,'/SVpipe.',sample,'.labeled.bed',sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}
label_PBNew<-function(sample, fam){
	pb_dat=read.table(paste(dat_path,'/PBNew.',sample,'.bed',sep=''))
	colnames(pb_dat)=c('CHR','POS','END','GT','SVTYPE')
	pb_dat[,ncol(pb_dat)+1]=pb_dat[,3]-pb_dat[,2]
	colnames(pb_dat)[ncol(pb_dat)]='SVLEN'
	pb_dat=order_bed(pb_dat)
	pb_dat=add_SV_Size(pb_dat)

        ill_cov=read.table(paste(s04_path, '/PB_New.',sample,'.bed.ILL_Cov', sep=''))
        ill_cov =  order_bed(ill_cov)
        pb_dat[,ncol(pb_dat)+1]=ill_cov[,7]
        colnames(pb_dat)[ncol(pb_dat)] = 'ILL_Cov'

	pb_dat$SVID=paste(pb_dat[,1],pb_dat[,2],pb_dat[,3],pb_dat[,4],sep='_')

	pb_dat=add_callset_supp(pb_dat, 'PBNew_GATK',sample, dat_path)
	pb_dat=add_callset_supp(pb_dat, 'PBNew_HGSV',sample, dat_path)
	pb_dat=add_callset_supp(pb_dat, 'PBNew_SVpipe',sample, dat_path)

	gc=read.table(paste(dat_path,'/genomic_context/PBNew.',sample,'.anno.bed',sep=''))
	gc=unique(order_bed(gc))
	gc[,ncol(gc)+1]=paste(gc[,1],gc[,2],gc[,3],gc[,4],sep='_')
	pb_dat[,ncol(pb_dat)+1]='-1'
	colnames(pb_dat)[ncol(pb_dat)]='GC'
	#pb_dat=pb_dat[,-9]
	pb_dat=unique(pb_dat)
	pb_dat[pb_dat$SVID %in% gc[,ncol(gc)],][,ncol(pb_dat)]=as.character(gc[,ncol(gc)-1])

	pb_dat$SVID=paste(pb_dat[,1],pb_dat[,2],pb_dat[,3],pb_dat[,5],sep='_')
	pb_dat[pb_dat$SVID %in% data.frame(table(pb_dat$SVID))[data.frame(table(pb_dat$SVID))[,2]>1,][,1] ,]$GT = 'HOM'
	pb_dat=unique(pb_dat)
	pb_dat=add_PB_Vali_to_PB(pb_dat,'PB_New',sample,'raw')
	pb_dat=add_PB_Vali_to_PB(pb_dat,'PB_New',sample,'h0')
	pb_dat=add_PB_Vali_to_PB(pb_dat,'PB_New',sample,'h1')

	write.table(pb_dat,paste(dat_path,'/PBNew.',sample,'.labeled.bed',sep=''), quote=F, sep='\t', col.names=T, row.names=F)
		
	}
sample_to_fam <-function(batch,sample){
	batch_dat=read.table(batch)
	out=batch_dat[batch_dat[,1]==sample,2]
	return(as.character(out))
	}


#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-c", "--callset"), type="character", default=NULL,
              help="name of callset", metavar="character"),
        make_option(c("-s", "--sample"), type="character", default=NULL,
              help="name of sample", metavar="character"),
       	make_option(c("-b", "--batch"), type="character", default=NULL,
              help="name of batch file", metavar="character"),
       	make_option(c("-d", "--data_path"), type="character", default='../data/raw_bed',
              help="path of input data", metavar="character"),
       	make_option(c("-p", "--pathof04"), type="character", default='../04_PB_Vali',
              help="name of batch file", metavar="character")
	);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

callset = opt$callset
sample = opt$sample
batch = opt$batch
dat_path=opt$data_path
s04_path=opt$pathof04

fam=sample_to_fam(batch, sample)
if(callset == 'GATK'){	
	label_GATK(sample, fam)
	}
if(callset == 'HGSV'){	
	label_HGSV(sample, fam)
	}
if(callset == 'SVpipe'){	
	label_SVpipe(sample, fam)
	}
if(callset == 'PBNew'){	
	label_PBNew(sample, fam)
	}


