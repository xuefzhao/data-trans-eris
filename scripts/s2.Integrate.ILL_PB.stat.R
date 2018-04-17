pb_vali_readin<-function(callset, sample){
	chs=read.table(paste(callset,sample,'PB_Supp.ILL_Supp.GC.bed',sep='.'), header=T)
	chs=chs[chs[,4]!='CPX',]
	vali_raw=read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.raw.bed.vapor',sep=''), header=T)
	vali_h0= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h0.bed.vapor',sep=''), header=T)
	vali_h1= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h1.bed.vapor',sep=''), header=T)

	vali_raw[,7]=as.character(vali_raw[,7])
	vali_raw[is.na(vali_raw$VaPoR_GT),][,7]=-1
	vali_h0[,7]=as.character(vali_h0[,7])
	vali_h0[is.na(vali_h0$VaPoR_GT),][,7]=-1
	vali_h1[,7]=as.character(vali_h1[,7])
	vali_h1[is.na(vali_h1$VaPoR_GT),][,7]=-1


	chs=chs[order(chs[,3]),]
	chs=chs[order(chs[,2]),]
	chs=chs[order(chs[,1]),]

	vali_raw=vali_raw[order(vali_raw[,3]),]
	vali_raw=vali_raw[order(vali_raw[,2]),]
	vali_raw=vali_raw[order(vali_raw[,1]),]

	vali_h0=vali_h0[order(vali_h0[,3]),]
	vali_h0=vali_h0[order(vali_h0[,2]),]
	vali_h0=vali_h0[order(vali_h0[,1]),]

	vali_h1=vali_h1[order(vali_h1[,3]),]
	vali_h1=vali_h1[order(vali_h1[,2]),]
	vali_h1=vali_h1[order(vali_h1[,1]),]

	chs[,ncol(chs)+1]=vali_raw$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h1'
	chs[,ncol(chs)+1]=vali_raw$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h1'
	return(chs)
	}
pb_HGSV_vali_readin<-function(callset, sample){
	chs=read.table(paste(callset,sample,'PB_Supp.ILL_Supp.GC.bed',sep='.'), header=T)
	chs=chs[chs[,4]!='CPX',]
	vali_raw=read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.raw.bed.vapor',sep=''), header=T)
	vali_h0= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h0.bed.vapor',sep=''), header=T)
	vali_h1= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h1.bed.vapor',sep=''), header=T)

	vali_raw[,7]=as.character(vali_raw[,7])
	vali_raw[is.na(vali_raw$VaPoR_GT),][,7]=-1
	vali_h0[,7]=as.character(vali_h0[,7])
	vali_h0[is.na(vali_h0$VaPoR_GT),][,7]=-1
	vali_h1[,7]=as.character(vali_h1[,7])
	vali_h1[is.na(vali_h1$VaPoR_GT),][,7]=-1

	chs[,ncol(chs)+1]=paste(chs[,1],chs[,2],chs[,3],sep='_')
	vali_raw[,ncol(vali_raw)+1]=paste(vali_raw[,1],vali_raw[,2],vali_raw[,3],sep='_')
	vali_h0[,ncol(vali_h0)+1]=paste(vali_h0[,1],vali_h0[,2],vali_h0[,3],sep='_')
	vali_h1[,ncol(vali_h1)+1]=paste(vali_h1[,1],vali_h1[,2],vali_h1[,3],sep='_')

	chs=chs[order(chs[,3]),]
	chs=chs[order(chs[,2]),]
	chs=chs[order(chs[,1]),]
	chs=chs[chs[,ncol(chs)] %in% vali_raw[,ncol(vali_raw)] ,]

	vali_raw=vali_raw[order(vali_raw[,3]),]
	vali_raw=vali_raw[order(vali_raw[,2]),]
	vali_raw=vali_raw[order(vali_raw[,1]),]

	vali_h0=vali_h0[order(vali_h0[,3]),]
	vali_h0=vali_h0[order(vali_h0[,2]),]
	vali_h0=vali_h0[order(vali_h0[,1]),]

	vali_h1=vali_h1[order(vali_h1[,3]),]
	vali_h1=vali_h1[order(vali_h1[,2]),]
	vali_h1=vali_h1[order(vali_h1[,1]),]

	chs[,ncol(chs)+1]=vali_raw$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h1'
	chs[,ncol(chs)+1]=vali_raw$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h1'
	return(chs)
	}
pb_PB_vali_readin<-function(callset, sample){
	chs=read.table(paste(callset,sample,'ILL_Supp.bed',sep='.'), header=T)
	chs=chs[chs[,4]!='CPX',]
	vali_raw=read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.raw.bed.vapor',sep=''), header=T)
	vali_h0= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h0.bed.vapor',sep=''), header=T)
	vali_h1= read.table(paste('../PB_Vali/',callset,'.',sample,'.PB_Vali.h1.bed.vapor',sep=''), header=T)

	vali_raw[,7]=as.character(vali_raw[,7])
	vali_raw[is.na(vali_raw$VaPoR_GT),][,7]=-1
	vali_h0[,7]=as.character(vali_h0[,7])
	vali_h0[is.na(vali_h0$VaPoR_GT),][,7]=-1
	vali_h1[,7]=as.character(vali_h1[,7])
	vali_h1[is.na(vali_h1$VaPoR_GT),][,7]=-1

	chs=chs[order(chs[,3]),]
	chs=chs[order(chs[,2]),]
	chs=chs[order(chs[,1]),]

	vali_raw=vali_raw[order(vali_raw[,3]),]
	vali_raw=vali_raw[order(vali_raw[,2]),]
	vali_raw=vali_raw[order(vali_raw[,1]),]

	vali_h0=vali_h0[order(vali_h0[,3]),]
	vali_h0=vali_h0[order(vali_h0[,2]),]
	vali_h0=vali_h0[order(vali_h0[,1]),]

	vali_h1=vali_h1[order(vali_h1[,3]),]
	vali_h1=vali_h1[order(vali_h1[,2]),]
	vali_h1=vali_h1[order(vali_h1[,1]),]

	chs[,ncol(chs)+1]=vali_raw$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GS
	colnames(chs)[ncol(chs)]='GS_h1'
	chs[,ncol(chs)+1]=vali_raw$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_raw'
	chs[,ncol(chs)+1]=vali_h0$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h0'
	chs[,ncol(chs)+1]=vali_h1$VaPoR_GT
	colnames(chs)[ncol(chs)]='GT_h1'

	ILL_cov=read.table(paste('PB_Uni.RO00.',sample,'.ILL_Cov.plus_anno',sep=''))
	ILL_cov[is.na(ILL_cov)]=-1
	chs[,ncol(chs)+1]=-1
	chs[chs$ID %in% ILL_cov[,7],][,ncol(chs)]=ILL_cov[,14]
	colnames(chs)[ncol(chs)]='ILL_Cov_median'

	return(chs)
	}
calcu_basic_stat_GATK<-function(gatk_chs, gatk_pur, gatk_yri, prefix){
	basic_stat=data.frame('HG00514'=0,'HG00733'=0,'NA19240'=0)
	row_num=1
	basic_stat[row_num,]=c(nrow(gatk_chs), nrow(gatk_pur), nrow(gatk_yri))
	rownames(basic_stat)[row_num]='NumSVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp=='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp=='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='PB_Supp_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp!='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp!='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp!='.5-1' & gatk_chs$HGSV_supp=='.5-1' & gatk_chs$SVpipe_supp!='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp!='.5-1' & gatk_pur$HGSV_supp=='.5-1' & gatk_pur$SVpipe_supp!='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp!='.5-1' & gatk_yri$HGSV_supp=='.5-1' & gatk_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.HGSV_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp!='.5-1' & gatk_chs$HGSV_supp!='.5-1' & gatk_chs$SVpipe_supp=='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp!='.5-1' & gatk_pur$HGSV_supp!='.5-1' & gatk_pur$SVpipe_supp=='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp!='.5-1' & gatk_yri$HGSV_supp!='.5-1' & gatk_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.SVpipe_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp!='.5-1' & gatk_chs$HGSV_supp=='.5-1' & gatk_chs$SVpipe_supp=='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp!='.5-1' & gatk_pur$HGSV_supp=='.5-1' & gatk_pur$SVpipe_supp=='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp!='.5-1' & gatk_yri$HGSV_supp=='.5-1' & gatk_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.HGSV_SVpipe_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(gatk_chs[gatk_chs$PBNew_supp!='.5-1' & gatk_chs$HGSV_supp!='.5-1' & gatk_chs$SVpipe_supp!='.5-1',]), 
							nrow(gatk_pur[gatk_pur$PBNew_supp!='.5-1' & gatk_pur$HGSV_supp!='.5-1' & gatk_pur$SVpipe_supp!='.5-1',]), 
							nrow(gatk_yri[gatk_yri$PBNew_supp!='.5-1' & gatk_yri$HGSV_supp!='.5-1' & gatk_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='GATK_Uni_SVs'

	basic_stat[,4]=prefix
	colnames(basic_stat)[4]='set'
	return(basic_stat)
	}
calcu_stat_GATK <-function(gatk_chs, gatk_pur, gatk_yri, prefix){
	basic_stat = calcu_basic_stat_GATK(gatk_chs, gatk_pur, gatk_yri, prefix)
	bi_parental = calcu_basic_stat_GATK(	gatk_chs[gatk_chs$inheri == 'fa_mo_pb',], 
											gatk_pur[gatk_pur$inheri == 'fa_mo_pb',], 
											gatk_yri[gatk_yri$inheri == 'fa_mo_pb',], prefix)
	uni_parental = calcu_basic_stat_GATK(   gatk_chs[gatk_chs$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_pur[gatk_pur$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_yri[gatk_yri$inheri%in%c('fa_pb','mo_pb'),], prefix)
	dn_parental = calcu_basic_stat_GATK(	gatk_chs[gatk_chs$inheri%in%c('denovo'),], 
											gatk_pur[gatk_pur$inheri%in%c('denovo'),], 
											gatk_yri[gatk_yri$inheri%in%c('denovo'),], prefix)
	basic_stat[,ncol(basic_stat)+1]='all'
	bi_parental[,ncol(bi_parental)+1]='bi_parental'
	uni_parental[,ncol(uni_parental)+1]='uni_parental'
	dn_parental[,ncol(dn_parental)+1]='denovo'
	out = rbind(basic_stat, bi_parental, uni_parental, dn_parental)
	colnames(out)[ncol(out)]='inheri'
	return(out)
	}
calcu_basic_stat_HGSV<-function(hgsv_chs, hgsv_pur, hgsv_yri, prefix){
	basic_stat=data.frame('HG00514'=0,'HG00733'=0,'NA19240'=0)
	row_num=1
	basic_stat[row_num,]=c(nrow(hgsv_chs), nrow(hgsv_pur), nrow(hgsv_yri))
	rownames(basic_stat)[row_num]='NumSVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp=='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp=='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='PB_Supp_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp!='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp!='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp!='.5-1' & hgsv_chs$GATK_supp=='.5-1' & hgsv_chs$SVpipe_supp!='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp!='.5-1' & hgsv_pur$GATK_supp=='.5-1' & hgsv_pur$SVpipe_supp!='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp!='.5-1' & hgsv_yri$GATK_supp=='.5-1' & hgsv_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.GATK_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp!='.5-1' & hgsv_chs$GATK_supp!='.5-1' & hgsv_chs$SVpipe_supp=='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp!='.5-1' & hgsv_pur$GATK_supp!='.5-1' & hgsv_pur$SVpipe_supp=='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp!='.5-1' & hgsv_yri$GATK_supp!='.5-1' & hgsv_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.SVpipe_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp!='.5-1' & hgsv_chs$GATK_supp=='.5-1' & hgsv_chs$SVpipe_supp=='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp!='.5-1' & hgsv_pur$GATK_supp=='.5-1' & hgsv_pur$SVpipe_supp=='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp!='.5-1' & hgsv_yri$GATK_supp=='.5-1' & hgsv_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.GATK_SVpipe_Supp'


	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(hgsv_chs[hgsv_chs$PBNew_supp!='.5-1' & hgsv_chs$GATK_supp!='.5-1' & hgsv_chs$SVpipe_supp!='.5-1',]), 
							nrow(hgsv_pur[hgsv_pur$PBNew_supp!='.5-1' & hgsv_pur$GATK_supp!='.5-1' & hgsv_pur$SVpipe_supp!='.5-1',]), 
							nrow(hgsv_yri[hgsv_yri$PBNew_supp!='.5-1' & hgsv_yri$GATK_supp!='.5-1' & hgsv_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='hgsv_Uni_SVs'

	basic_stat[,4]=prefix
	colnames(basic_stat)[4]='set'
	return(basic_stat)
	}
calcu_stat_HGSV <-function(gatk_chs, gatk_pur, gatk_yri, prefix){
	basic_stat = calcu_basic_stat_HGSV(gatk_chs, gatk_pur, gatk_yri, prefix)
	bi_parental = calcu_basic_stat_HGSV(	gatk_chs[gatk_chs$inheri == 'fa_mo_pb',], 
											gatk_pur[gatk_pur$inheri == 'fa_mo_pb',], 
											gatk_yri[gatk_yri$inheri == 'fa_mo_pb',], prefix)
	uni_parental = calcu_basic_stat_HGSV(   gatk_chs[gatk_chs$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_pur[gatk_pur$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_yri[gatk_yri$inheri%in%c('fa_pb','mo_pb'),], prefix)
	dn_parental = calcu_basic_stat_HGSV(	gatk_chs[gatk_chs$inheri%in%c('denovo'),], 
											gatk_pur[gatk_pur$inheri%in%c('denovo'),], 
											gatk_yri[gatk_yri$inheri%in%c('denovo'),], prefix)
	basic_stat[,ncol(basic_stat)+1]='all'
	bi_parental[,ncol(bi_parental)+1]='bi_parental'
	uni_parental[,ncol(uni_parental)+1]='uni_parental'
	dn_parental[,ncol(dn_parental)+1]='denovo'
	out = rbind(basic_stat, bi_parental, uni_parental, dn_parental)
	colnames(out)[ncol(out)]='inheri'
	return(out)
	}
calcu_basic_stat_SVpipe<-function(svpipe_chs, svpipe_pur, svpipe_yri, prefix){
	basic_stat=data.frame('HG00514'=0,'HG00733'=0,'NA19240'=0)
	row_num=1
	basic_stat[row_num,]=c(nrow(svpipe_chs), nrow(svpipe_pur), nrow(svpipe_yri))
	rownames(basic_stat)[row_num]='NumSVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(svpipe_chs[svpipe_chs$PBNew_supp=='.5-1',]), 
							nrow(svpipe_pur[svpipe_pur$PBNew_supp=='.5-1',]), 
							nrow(svpipe_yri[svpipe_yri$PBNew_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='PB_Supp_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(svpipe_chs[svpipe_chs$PBNew_supp!='.5-1',]), 
							nrow(svpipe_pur[svpipe_pur$PBNew_supp!='.5-1',]), 
							nrow(svpipe_yri[svpipe_yri$PBNew_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(svpipe_chs[svpipe_chs$PBNew_supp!='.5-1' & svpipe_chs$GATK_supp=='.5-1' & svpipe_chs$HGSV_supp!='.5-1',]), 
							nrow(svpipe_pur[svpipe_pur$PBNew_supp!='.5-1' & svpipe_pur$GATK_supp=='.5-1' & svpipe_pur$HGSV_supp!='.5-1',]), 
							nrow(svpipe_yri[svpipe_yri$PBNew_supp!='.5-1' & svpipe_yri$GATK_supp=='.5-1' & svpipe_yri$HGSV_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.GATK_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(svpipe_chs[svpipe_chs$PBNew_supp!='.5-1' & svpipe_chs$GATK_supp!='.5-1' & svpipe_chs$HGSV_supp=='.5-1',]), 
							nrow(svpipe_pur[svpipe_pur$PBNew_supp!='.5-1' & svpipe_pur$GATK_supp!='.5-1' & svpipe_pur$HGSV_supp=='.5-1',]), 
							nrow(svpipe_yri[svpipe_yri$PBNew_supp!='.5-1' & svpipe_yri$GATK_supp!='.5-1' & svpipe_yri$HGSV_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.HGSV_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c( nrow(svpipe_chs[svpipe_chs$PBNew_supp!='.5-1' & svpipe_chs$GATK_supp=='.5-1' & svpipe_chs$HGSV_supp=='.5-1',]), 
							nrow(svpipe_pur[svpipe_pur$PBNew_supp!='.5-1' & svpipe_pur$GATK_supp=='.5-1' & svpipe_pur$HGSV_supp=='.5-1',]), 
							nrow(svpipe_yri[svpipe_yri$PBNew_supp!='.5-1' & svpipe_yri$GATK_supp=='.5-1' & svpipe_yri$HGSV_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Uni_SVs.GATK_HGSV_Supp'

	row_num=row_num+1
	basic_stat[row_num,]=c(nrow(svpipe_chs[svpipe_chs$PBNew_supp!='.5-1' & svpipe_chs$GATK_supp!='.5-1' & svpipe_chs$HGSV_supp!='.5-1',]), nrow(svpipe_pur[svpipe_pur$PBNew_supp!='.5-1' & svpipe_pur$GATK_supp!='.5-1' & svpipe_pur$HGSV_supp!='.5-1',]), nrow(svpipe_yri[svpipe_yri$PBNew_supp!='.5-1' & svpipe_yri$GATK_supp!='.5-1' & svpipe_yri$HGSV_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='svpipe_Uni_SVs'

	basic_stat[,4]=prefix
	colnames(basic_stat)[4]='set'
	return(basic_stat)
	}
calcu_stat_SVpipe <-function(gatk_chs, gatk_pur, gatk_yri, prefix){
	basic_stat = calcu_basic_stat_SVpipe(gatk_chs, gatk_pur, gatk_yri, prefix)
	bi_parental = calcu_basic_stat_SVpipe(	gatk_chs[gatk_chs$inheri == 'fa_mo_pb',], 
											gatk_pur[gatk_pur$inheri == 'fa_mo_pb',], 
											gatk_yri[gatk_yri$inheri == 'fa_mo_pb',], prefix)
	uni_parental = calcu_basic_stat_SVpipe(   gatk_chs[gatk_chs$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_pur[gatk_pur$inheri%in%c('fa_pb','mo_pb'),], 
											gatk_yri[gatk_yri$inheri%in%c('fa_pb','mo_pb'),], prefix)
	dn_parental = calcu_basic_stat_SVpipe(	gatk_chs[gatk_chs$inheri%in%c('denovo'),], 
											gatk_pur[gatk_pur$inheri%in%c('denovo'),], 
											gatk_yri[gatk_yri$inheri%in%c('denovo'),], prefix)
	basic_stat[,ncol(basic_stat)+1]='all'
	bi_parental[,ncol(bi_parental)+1]='bi_parental'
	uni_parental[,ncol(uni_parental)+1]='uni_parental'
	dn_parental[,ncol(dn_parental)+1]='denovo'
	out = rbind(basic_stat, bi_parental, uni_parental, dn_parental)
	colnames(out)[ncol(out)]='inheri'
	return(out)
	}
calcu_basic_stat_PB<-function(pb_chs, pb_pur, pb_yri,prefix){
	basic_stat=data.frame('HG00514'=0,'HG00733'=0,'NA19240'=0)
	row_num=1
	basic_stat[row_num,]=c(nrow(pb_chs), nrow(pb_pur), nrow(pb_yri))
	rownames(basic_stat)[row_num]='NumSVs'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp!='.5-1' & pb_chs$HGSV_supp!='.5-1' & pb_chs$SVpipe_supp!='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp!='.5-1' & pb_pur$HGSV_supp!='.5-1' & pb_pur$SVpipe_supp!='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp!='.5-1' & pb_yri$HGSV_supp!='.5-1' & pb_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='PB_Uni'

	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp=='.5-1' & pb_chs$HGSV_supp!='.5-1' & pb_chs$SVpipe_supp!='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp=='.5-1' & pb_pur$HGSV_supp!='.5-1' & pb_pur$SVpipe_supp!='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp=='.5-1' & pb_yri$HGSV_supp!='.5-1' & pb_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_GATK'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp!='.5-1' & pb_chs$HGSV_supp=='.5-1' & pb_chs$SVpipe_supp!='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp!='.5-1' & pb_pur$HGSV_supp=='.5-1' & pb_pur$SVpipe_supp!='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp!='.5-1' & pb_yri$HGSV_supp=='.5-1' & pb_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_HGSV'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp!='.5-1' & pb_chs$HGSV_supp!='.5-1' & pb_chs$SVpipe_supp=='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp!='.5-1' & pb_pur$HGSV_supp!='.5-1' & pb_pur$SVpipe_supp=='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp!='.5-1' & pb_yri$HGSV_supp!='.5-1' & pb_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_SVpipe'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp=='.5-1' & pb_chs$HGSV_supp=='.5-1' & pb_chs$SVpipe_supp!='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp=='.5-1' & pb_pur$HGSV_supp=='.5-1' & pb_pur$SVpipe_supp!='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp=='.5-1' & pb_yri$HGSV_supp=='.5-1' & pb_yri$SVpipe_supp!='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_GATK_HGSV'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp=='.5-1' & pb_chs$HGSV_supp!='.5-1' & pb_chs$SVpipe_supp=='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp=='.5-1' & pb_pur$HGSV_supp!='.5-1' & pb_pur$SVpipe_supp=='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp=='.5-1' & pb_yri$HGSV_supp!='.5-1' & pb_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_GATK_SVpipe'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp!='.5-1' & pb_chs$HGSV_supp=='.5-1' & pb_chs$SVpipe_supp=='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp!='.5-1' & pb_pur$HGSV_supp=='.5-1' & pb_pur$SVpipe_supp=='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp!='.5-1' & pb_yri$HGSV_supp=='.5-1' & pb_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_HGSV_SVpipe'
	row_num=1+row_num
	basic_stat[row_num,]=c( nrow(pb_chs[pb_chs$GATK_supp=='.5-1' & pb_chs$HGSV_supp=='.5-1' & pb_chs$SVpipe_supp=='.5-1',]), 
							nrow(pb_pur[pb_pur$GATK_supp=='.5-1' & pb_pur$HGSV_supp=='.5-1' & pb_pur$SVpipe_supp=='.5-1',]), 
							nrow(pb_yri[pb_yri$GATK_supp=='.5-1' & pb_yri$HGSV_supp=='.5-1' & pb_yri$SVpipe_supp=='.5-1',]))
	rownames(basic_stat)[row_num]='ILL_Supp_GATK_HGSV_SVpipe'
	basic_stat[,4]=prefix
	return(basic_stat)
	}
calcu_stat_PB<-function(pb_chs, pb_pur, pb_yri,prefix){
	basic_stat = calcu_basic_stat_PB(pb_chs, pb_pur, pb_yri,prefix)
	hom_del = calcu_basic_stat_PB(pb_chs[pb_chs$GT=='HOM' & pb_chs$SVTYPE=='DEL' ,], pb_pur[pb_pur$GT=='HOM' & pb_pur$SVTYPE=='DEL' ,],pb_yri[pb_yri$GT=='HOM' & pb_yri$SVTYPE=='DEL' ,],prefix)
	het_del = calcu_basic_stat_PB(pb_chs[pb_chs$GT!='HOM' & pb_chs$SVTYPE=='DEL' ,], pb_pur[pb_pur$GT!='HOM' & pb_pur$SVTYPE=='DEL' ,],pb_yri[pb_yri$GT!='HOM' & pb_yri$SVTYPE=='DEL' ,],prefix)
	hom_ins = calcu_basic_stat_PB(pb_chs[pb_chs$GT=='HOM' & pb_chs$SVTYPE=='INS' ,], pb_pur[pb_pur$GT=='HOM' & pb_pur$SVTYPE=='INS' ,],pb_yri[pb_yri$GT=='HOM' & pb_yri$SVTYPE=='INS' ,],prefix)
	het_ins = calcu_basic_stat_PB(pb_chs[pb_chs$GT!='HOM' & pb_chs$SVTYPE=='INS' ,], pb_pur[pb_pur$GT!='HOM' & pb_pur$SVTYPE=='INS' ,],pb_yri[pb_yri$GT!='HOM' & pb_yri$SVTYPE=='INS' ,],prefix)
	basic_stat[,ncol(basic_stat)+1]='all'
	hom_del[,ncol(hom_del)+1]='hom_del'
	het_del[,ncol(het_del)+1]='het_del'
	hom_ins[,ncol(hom_ins)+1]='hom_ins'
	het_ins[,ncol(het_ins)+1]='het_ins'
	out=rbind(basic_stat, hom_del, het_del, hom_ins, het_ins)
	colnames(out)[ncol(out)-1]='set'
	colnames(out)[ncol(out)]='gt_svtype'
	return(out)
	}

reformat_stat<-function(sample_name,callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat){
	if (sample_name=='HG00514'){ col_Num=1 }
	if (sample_name=='HG00733'){ col_Num=2 }
	if (sample_name=='NA19240'){ col_Num=3 }
       GATK_chs_stat=cbind(basic_stat[,col_Num], PB_raw_vali_stat[,col_Num], PB_assm_vali_stat[,col_Num], PB_raw_only_vali_stat[,col_Num], PB_assm_only_vali_stat[,col_Num], PB_raw_assm_vali_stat[,col_Num], PB_no_vali_stat[,col_Num],basic_stat[,ncol(basic_stat)],sample_name)
        colnames(GATK_chs_stat)=c(callset,'PB_raw_Vali','PB_assem_Vali','PB_raw_only_Vali', 'PB_assm_only_Vali', 'PB_raw_assm_Vali', 'PB_not_Vali', 'inheri','sample')
        rownames(GATK_chs_stat)=rownames(basic_stat)
        return(GATK_chs_stat)
        }

labeled_data_readin<-function(callset,sample_name, dat_path){
	dat_chs=read.table(paste(dat_path, '/',callset,'.',sample_name,'.labeled.bed',sep=''), header=T)
	dat_chs[is.na(dat_chs)]=-1
	return(dat_chs)
	}

GATK_dat_process<-function(dat_path){
	callset	= 'GATK'
	dat_chs = labeled_data_readin(callset, 'HG00514', dat_path)
	dat_pur = labeled_data_readin(callset, 'HG00733', dat_path)
	dat_yri = labeled_data_readin(callset, 'NA19240', dat_path)

	basic_stat = calcu_stat_GATK(dat_chs, dat_pur, dat_yri, paste(callset, 'all', sep='.'))
	PB_raw_vali_stat = calcu_stat_GATK(	dat_chs[dat_chs$PB_vali_raw %in% c('0/1','1/1'),], 	dat_pur[dat_pur$PB_vali_raw %in% c('0/1','1/1'),], 	dat_yri[dat_yri$PB_vali_raw %in% c('0/1','1/1'),], paste(callset,'PB_raw_Vali',sep='.'))
	PB_assm_vali_stat = calcu_stat_GATK(	unique(rbind(dat_chs[dat_chs$PB_vali_h0 %in% c('0/1','1/1'), ],dat_chs[dat_chs$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_pur[dat_pur$PB_vali_h0 %in% c('0/1','1/1'), ],dat_pur[dat_pur$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_yri[dat_yri$PB_vali_h0 %in% c('0/1','1/1'), ],dat_yri[dat_yri$PB_vali_h1 %in% c('0/1','1/1'), ])), paste(callset, 'PB_assem_Vali',sep='.'))
	PB_raw_only_vali_stat = calcu_stat_GATK(	dat_chs[ dat_chs$PB_vali_raw %in% c('0/1','1/1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[ dat_pur$PB_vali_raw %in% c('0/1','1/1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[ dat_yri$PB_vali_raw %in% c('0/1','1/1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_raw_only_Vali',sep='.'))
	PB_assm_only_vali_stat = calcu_stat_GATK(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/0', '-1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/0', '-1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/0', '-1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_assm_only_Vali',sep='.'))
	PB_raw_assm_vali_stat = calcu_stat_GATK(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/1','1/1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/1','1/1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/1','1/1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_raw_assm_Vali',sep='.'))
	PB_no_vali_stat = calcu_stat_GATK(	dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_not_Vali',sep='.'))

	dat_stats=rbind(basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat)
	dat_stats_2=rbind(	reformat_stat('HG00514',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('HG00733',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('NA19240',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat))
	dat_stats_2=cbind(rownames(dat_stats_2), dat_stats_2)
	colnames(dat_stats_2)[1]='names'
	write.table(dat_stats,   paste('stats/ILL_vs_PB.',callset,'_V1.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=T)
	write.table(dat_stats_2, paste('stats/ILL_vs_PB.',callset,'_V2.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}

HGSV_dat_process<-function(dat_path){
	callset	= 'HGSV'
	dat_chs = labeled_data_readin(callset, 'HG00514', dat_path)
	dat_pur = labeled_data_readin(callset, 'HG00733', dat_path)
	dat_yri = labeled_data_readin(callset, 'NA19240', dat_path)

	basic_stat = calcu_stat_HGSV(dat_chs, dat_pur, dat_yri, paste(callset, 'all', sep='.'))
	PB_raw_vali_stat = calcu_stat_HGSV(	dat_chs[dat_chs$PB_vali_raw %in% c('0/1','1/1'),], 	dat_pur[dat_pur$PB_vali_raw %in% c('0/1','1/1'),], 	dat_yri[dat_yri$PB_vali_raw %in% c('0/1','1/1'),], paste(callset,'PB_raw_Vali',sep='.'))
	PB_assm_vali_stat = calcu_stat_HGSV(	unique(rbind(dat_chs[dat_chs$PB_vali_h0 %in% c('0/1','1/1'), ],dat_chs[dat_chs$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_pur[dat_pur$PB_vali_h0 %in% c('0/1','1/1'), ],dat_pur[dat_pur$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_yri[dat_yri$PB_vali_h0 %in% c('0/1','1/1'), ],dat_yri[dat_yri$PB_vali_h1 %in% c('0/1','1/1'), ])), paste(callset, 'PB_assem_Vali',sep='.'))
	PB_raw_only_vali_stat = calcu_stat_HGSV(	dat_chs[ dat_chs$PB_vali_raw %in% c('0/1','1/1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[ dat_pur$PB_vali_raw %in% c('0/1','1/1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[ dat_yri$PB_vali_raw %in% c('0/1','1/1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_raw_only_Vali',sep='.'))
	PB_assm_only_vali_stat = calcu_stat_HGSV(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/0', '-1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/0', '-1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/0', '-1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_assm_only_Vali',sep='.'))
	PB_raw_assm_vali_stat = calcu_stat_HGSV(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/1','1/1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/1','1/1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/1','1/1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_raw_assm_Vali',sep='.'))
	PB_no_vali_stat = calcu_stat_HGSV(	dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_not_Vali',sep='.'))

	dat_stats=rbind(basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat)
	dat_stats_2=rbind(	reformat_stat('HG00514',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('HG00733',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('NA19240',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat))
	dat_stats_2=cbind(rownames(dat_stats_2), dat_stats_2)
	colnames(dat_stats_2)[1]='names'
	write.table(dat_stats,   paste('stats/ILL_vs_PB.',callset,'_V1.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=T)
	write.table(dat_stats_2, paste('stats/ILL_vs_PB.',callset,'_V2.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}

SVpipe_dat_process<-function(dat_path){
	callset = 'SVpipe'
	dat_chs = labeled_data_readin(callset, 'HG00514', dat_path)
	dat_pur = labeled_data_readin(callset, 'HG00733', dat_path)
	dat_yri = labeled_data_readin(callset, 'NA19240', dat_path)

	basic_stat = calcu_stat_SVpipe(dat_chs, dat_pur, dat_yri, paste(callset, 'all', sep='.'))
	PB_raw_vali_stat = calcu_stat_SVpipe(	dat_chs[dat_chs$PB_vali_raw %in% c('0/1','1/1'),], 	dat_pur[dat_pur$PB_vali_raw %in% c('0/1','1/1'),], 	dat_yri[dat_yri$PB_vali_raw %in% c('0/1','1/1'),], paste(callset,'PB_raw_Vali',sep='.'))
	PB_assm_vali_stat = calcu_stat_SVpipe(	unique(rbind(dat_chs[dat_chs$PB_vali_h0 %in% c('0/1','1/1'), ],dat_chs[dat_chs$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_pur[dat_pur$PB_vali_h0 %in% c('0/1','1/1'), ],dat_pur[dat_pur$PB_vali_h1 %in% c('0/1','1/1'), ])),	unique(rbind(dat_yri[dat_yri$PB_vali_h0 %in% c('0/1','1/1'), ],dat_yri[dat_yri$PB_vali_h1 %in% c('0/1','1/1'), ])), paste(callset, 'PB_assem_Vali',sep='.'))
	PB_raw_only_vali_stat = calcu_stat_SVpipe(	dat_chs[ dat_chs$PB_vali_raw %in% c('0/1','1/1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[ dat_pur$PB_vali_raw %in% c('0/1','1/1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[ dat_yri$PB_vali_raw %in% c('0/1','1/1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_raw_only_Vali',sep='.'))
	PB_assm_only_vali_stat = calcu_stat_SVpipe(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/0', '-1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/0', '-1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/0', '-1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_assm_only_Vali',sep='.'))
	PB_raw_assm_vali_stat = calcu_stat_SVpipe(	unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/1','1/1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/1','1/1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])),	unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/1','1/1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_raw_assm_Vali',sep='.'))
	PB_no_vali_stat = calcu_stat_SVpipe(	dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),],	dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),],	dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_not_Vali',sep='.'))

	dat_stats=rbind(basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat)
	dat_stats_2=rbind(	reformat_stat('HG00514',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('HG00733',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('NA19240',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat))
	dat_stats_2=cbind(rownames(dat_stats_2), dat_stats_2)
	colnames(dat_stats_2)[1]='names'
	write.table(dat_stats,   paste('stats/ILL_vs_PB.',callset,'_V1.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=T)
	write.table(dat_stats_2, paste('stats/ILL_vs_PB.',callset,'_V2.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}

PBNew_dat_process<-function(dat_path){
	callset = 'PBNew'
	dat_chs = labeled_data_readin(callset, 'HG00514', dat_path)
	dat_pur = labeled_data_readin(callset, 'HG00733', dat_path)
	dat_yri = labeled_data_readin(callset, 'NA19240', dat_path)

	basic_stat = calcu_stat_PB(dat_chs, dat_pur, dat_yri, paste(callset, 'all', sep='.'))
	PB_raw_vali_stat = calcu_stat_PB( dat_chs[dat_chs$PB_vali_raw %in% c('0/1','1/1'),],  dat_pur[dat_pur$PB_vali_raw %in% c('0/1','1/1'),],  dat_yri[dat_yri$PB_vali_raw %in% c('0/1','1/1'),], paste(callset,'PB_raw_Vali',sep='.'))
	PB_assm_vali_stat = calcu_stat_PB( unique(rbind(dat_chs[dat_chs$PB_vali_h0 %in% c('0/1','1/1'), ],dat_chs[dat_chs$PB_vali_h1 %in% c('0/1','1/1'), ])), unique(rbind(dat_pur[dat_pur$PB_vali_h0 %in% c('0/1','1/1'), ],dat_pur[dat_pur$PB_vali_h1 %in% c('0/1','1/1'), ])), unique(rbind(dat_yri[dat_yri$PB_vali_h0 %in% c('0/1','1/1'), ],dat_yri[dat_yri$PB_vali_h1 %in% c('0/1','1/1'), ])), paste(callset, 'PB_assem_Vali',sep='.'))
	PB_raw_only_vali_stat = calcu_stat_PB( dat_chs[ dat_chs$PB_vali_raw %in% c('0/1','1/1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),], dat_pur[ dat_pur$PB_vali_raw %in% c('0/1','1/1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),], dat_yri[ dat_yri$PB_vali_raw %in% c('0/1','1/1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_raw_only_Vali',sep='.'))
	PB_assm_only_vali_stat = calcu_stat_PB( unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/0', '-1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])), unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/0', '-1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])), unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/0', '-1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_assm_only_Vali',sep='.'))
	PB_raw_assm_vali_stat = calcu_stat_PB( unique(rbind(dat_chs[dat_chs$PB_vali_raw %in%c('0/1','1/1') & dat_chs$PB_vali_h0%in% c('0/1','1/1'),],dat_chs[dat_chs$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_chs$PB_vali_h1) & dat_chs$PB_vali_h1%in% c('0/1','1/1'),])), unique(rbind(dat_pur[dat_pur$PB_vali_raw %in%c('0/1','1/1') & dat_pur$PB_vali_h0%in% c('0/1','1/1'),],dat_pur[dat_pur$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_pur$PB_vali_h1) & dat_pur$PB_vali_h1%in% c('0/1','1/1'),])), unique(rbind(dat_yri[dat_yri$PB_vali_raw %in%c('0/1','1/1') & dat_yri$PB_vali_h0%in% c('0/1','1/1'),],dat_yri[dat_yri$PB_vali_raw  %in%c('0/1','1/1') & !is.na(dat_yri$PB_vali_h1) & dat_yri$PB_vali_h1%in% c('0/1','1/1'),])), paste(callset, 'PB_raw_assm_Vali',sep='.'))
	PB_no_vali_stat = calcu_stat_PB( dat_chs[dat_chs$PB_vali_raw  %in%c('0/0', '-1') & dat_chs$PB_vali_h0 %in%c('0/0', '-1') & dat_chs$PB_vali_h1 %in%c('0/0', '-1'),], dat_pur[dat_pur$PB_vali_raw  %in%c('0/0', '-1') & dat_pur$PB_vali_h0 %in%c('0/0', '-1') & dat_pur$PB_vali_h1 %in%c('0/0', '-1'),], dat_yri[dat_yri$PB_vali_raw  %in%c('0/0', '-1') & dat_yri$PB_vali_h0 %in%c('0/0', '-1') & dat_yri$PB_vali_h1 %in%c('0/0', '-1'),], paste(callset, 'PB_not_Vali',sep='.'))

	dat_stats=rbind(basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat)
	dat_stats_2=rbind(	reformat_stat('HG00514',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('HG00733',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat), 
						reformat_stat('NA19240',callset,basic_stat, PB_raw_vali_stat, PB_assm_vali_stat, PB_raw_only_vali_stat, PB_assm_only_vali_stat, PB_raw_assm_vali_stat, PB_no_vali_stat))
	dat_stats_2=cbind(rownames(dat_stats_2), dat_stats_2)
	colnames(dat_stats_2)[1]='names'
	write.table(dat_stats,   paste('stats/ILL_vs_PB.',callset,'_V1.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=T)
	write.table(dat_stats_2, paste('stats/ILL_vs_PB.',callset,'_V2.stat' , sep=''), quote=F, sep='\t', col.names=T, row.names=F)
	}



#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-c", "--callset"), type="character", default=NULL,
              help="name of callset", metavar="character"),
        make_option(c("-d", "--data_path"), type="character", default='../data/raw_bed',
              help="path of input data", metavar="character")
        );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

callset = opt$callset
dat_path=opt$data_path
if(callset=='HGSV'){
	HGSV_dat_process(dat_path)
	}

if(callset=='GATK'){
	GATK_dat_process(dat_path)
	}

if(callset=='SVpipe'){
	SVpipe_dat_process(dat_path)
	}

if(callset=='PBNew'){
	PBNew_dat_process(dat_path)
	}



