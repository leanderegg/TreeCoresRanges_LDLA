#June 28, 2012
#Goal: Examine growth-climate and growth-competition relationships in sapling and adult trees at Mt. Rainier to answer the question: How do climate and competition interact to determine tree range limits?
#To answer this question,fit models of:
#1) Climate-growth relationships across elevation in saplings & Adults
	#growth measure: 
	####a. basal area increment (BAI)
	####b. relative growth rate (BAI/total BA)
#2) Competition-growth relationships across elevation in saplings & Adults
#Use coefficients or average growth across elevation to show graph as a summary.

#Directory:"/Users/ailene_kane/Documents/UW/Research/Mount Rainier/saplings/Summer2012Analyses"
#read in files for adult trees
treegrowth<-read.csv("Ailenecode/SouthSidecores2008-2011.csv", header=TRUE)
treegrowth[1:5,1:10]
dim(treegrowth)
unique(treegrowth[,4])#612 cored trees
speciesfoc<-unique(treegrowth$Species)#6 species
standfoc<-unique(treegrowth$Stand)#9 stands
rawdat<-treegrowth[,c(1:4,6:7,10:170)] #only columns we need
rawdat[1:5,1:10]
tags<-unique(rawdat[,4])
inddat<-matrix(NA, nrow=length(tags), ncol=6)
dimnames(inddat)<-list(c(), c("stand","elev","species","tag","DBH","lastyear"))
rngdat<-c()

for(i in 1:length(tags)){
	tmpdat<-rawdat[rawdat[,4]==tags[i],] #loop through trees
	lastring<-min(tmpdat[,6])-1 # drop the ring from the year collected?
	if(lastring==2007){firstyr<-9}
	if(lastring==2008){firstyr<-8}
	lastyr<-c()
	for(j in 1:dim(tmpdat)[1]){ # loop through cores to find the last year of each
		rngs<-tmpdat[j,7:dim(tmpdat)[2]] # get # of years
		tmp<-which(is.na(rngs)==TRUE) # get years beyond what's measured
		tmp<-tmp[tmp>3] #make sure there are at least 3 years removed
		if(length(tmp)==0){mintmp<-161}
		if(length(tmp)>0){mintmp<-min(tmp)}
		lastyr<-c(lastyr,mintmp+4)
	}
	lastyr<-min(lastyr) # find the shortest core 
  
	avgrng<-colMeans(tmpdat[,firstyr:lastyr]) # average cores
  
	for(j in 1:5){
		inddat[i,j]<-as.character(tmpdat[1,j])}
	inddat[i,6]<-lastring
	if(lastring==2007){
		treerng<-c(NA,avgrng,rep(NA, times=dim(tmpdat)[2]-lastyr))} # putting the NAs back into the data?
	if(lastring==2008){
		treerng<-c(avgrng,rep(NA, times=dim(tmpdat)[2]-lastyr))}
	rngdat<-rbind(rngdat,treerng)
}

# Whew. I think all that was significantly harder to follow than my code- LDLA

dimnames(rngdat)<-list(c(), seq(2008,1849))
dat<-cbind(inddat,rngdat[,1:159]) #merge data together
dat.df<-as.data.frame(cbind(inddat,rngdat[,1:159]))
rownames(rngdat)<-dat[,3]
dat.df[1:5,1:10]#this data frame contains info on species, stand, dbh (in cm), annual ring widths (in mm) going back to 1850

#convert annual ring widths to BAI and RGR:
inddat<-c()
baidat<-c()
rgrdat<-c()
for(i in 1:length(tags)){
	tmptreedat<-dat[dat[,4]==tags[i],]
	lastring<-as.numeric(tmptreedat[6])-1#skip last year of growth- not accurate (sampled part way through growing season)
	if(lastring==2006){lastyr<-9}
	if(lastring==2007){lastyr<-8}
	firstyr<-min(which(is.na(tmptreedat)))
	years<-length(lastyr:firstyr)-1#skip earliest year of growth, as may be inaccurate
	indbai<-c()
	indrgr<-c()
	for (j in 1:years){
		t1<-j-1
		t2<-j
		radius1<-((as.numeric(tmptreedat[5])*10)/2)-(sum(as.numeric(tmptreedat[(lastyr-1):(lastyr-1+t1)])))#radius in mm at time t1 (taking into account that the year core was taken is not a full ring width substract the growth from the year of collection plus growth in subsequent years from dbh to get this #
		#print(radius1)
		radius2<-radius1-as.numeric(tmptreedat[(lastyr-1+t2)])#radius at time t2, in mm
		#print(radius2)
		ba1<-(pi*radius1*radius1)
		#print(ba1)
		ba2<-(pi*radius2*radius2)
		#print(ba2)
		bai<-ba1-ba2
		rgr<-(ba1-ba2)/ba1#relative growth rate=bai/ba
		#print(bai)
		indbai<-c(indbai,bai)#tacks on next) oldest year of bai to end of indbai
		indrgr<-c(indrgr,rgr)#tacks on next) oldest year of rgr to end of indrgr
		}
		#bai data
	if(lastring==2006){indbai2<-c(NA,indbai)}
	if(lastring==2007){indbai2<-indbai}
		indbai3<-c(indbai2,rep(NA,times=(158-(length(indbai2)))))
		baidat<-rbind(baidat,indbai3)
		indtreedat<-tmptreedat[1:6]
		inddat<-rbind(inddat,indtreedat)
		#rgr data
		if(lastring==2006){indrgr2<-c(NA,indrgr)}
	if(lastring==2007){indrgr2<-indrgr}
		indrgr3<-c(indrgr2,rep(NA,times=(158-(length(indrgr2)))))
		rgrdat<-rbind(rgrdat,indrgr3)
	}
		colnames(baidat)<-seq(from=2007, to=1850)
		colnames(rgrdat)<-seq(from=2007, to=1850)
		treebaidat<-cbind(inddat,baidat)
		treergrdat<-cbind(inddat,rgrdat)
		
treebaidat[1:100,1:8]
treergrdat[1:5,1:17]
dat[1:5,1:10]
	plot(dat[150,8:58],treebaidat[150,8:58])

#read in files for saplings
sapdat<-read.csv("sapdat2012.csv",header=TRUE)#includes annual growth (avg ring width (mm)) for ABAM and TSME saplings at and some stands. Need to add gap status. 
dim(sapdat)
sapdat<-sapdat[1:1158,]#no growth data on GS42, the last sapling listed there...not sure why!
sapdat<-sapdat[sapdat$TreeName!="ES34",]#gap status not known for this individual. Remove from dataset
sapdat<-sapdat[sapdat$TreeName!="ES22",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="ES36",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="AS41",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="FS01",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="FS36",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="FS39",]#fewer than 4 paths, so no growth data for this tree 
sapdat<-sapdat[sapdat$TreeName!="GS20",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="GS21",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="HS06",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="HS35",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="HS49",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS06",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS07",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS18",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS30",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS35",]#fewer than 4 paths, so no growth data for this tree
sapdat<-sapdat[sapdat$TreeName!="PS39",]#fewer than 4 paths, so no growth data for this tree


dim(sapdat)
sapdat<-sapdat[,1:120]#just go back to 1910

#Get Basal Area Increment and relative growth rate (RGR=BAI/BA)from ring widths of 4 paths. This code takes a while to run.

stands<-unique(sapdat[,2])
stands
sapnum<-unique(sapdat[,4])#new sap number
sapnum#272 different saplings
gapstatus<-unique(sapdat[,10])
gapstatus#gap, not gap
species<-unique(sapdat[,5])
species#Tsme, Abam, Tshe
inddat<-c()
baidat<-c()
rgrdat<-c()
for(i in 1:length(sapnum)){
	tmpdat<-sapdat[sapdat[,4]==sapnum[i],]
	tmpdat2<-cbind(tmpdat[,2],tmpdat[,4:5],tmpdat[,9:10],tmpdat[,15])
	lastring<-min(tmpdat[,15])-1#skip last year of growth- not accurate (sampled part way through growing season)
	if(lastring==2007){lastyr<-23}
	if(lastring==2008){lastyr<-22}
	if(lastring==2009){lastyr<-21}
	if(lastring==2010){lastyr<-20}
	firstyr<-lastyr+(tmpdat[1,16]-2)
	if(firstyr>120){firstyr<-120}
	lastyr<-min(lastyr)
	ringcount<-tmpdat[1,16]
	years<-length(lastyr:firstyr)-2
	indbai<-c()
	indrgr<-c()
	for (j in 1:years){
		t1<-j-1
		t2<-j
		radii1<-rowSums(tmpdat[,(lastyr+t1):firstyr], na.rm=TRUE)#radii at time t1
		#print(radii1)
		radii2<-rowSums(tmpdat[,(lastyr+t2):firstyr], na.rm=TRUE)#radii at time t2(which equals t-1)
		#print(radii2)
		ba1<-(pi*radii1[1]*radii1[2]/4)+(pi*radii1[2]*radii1[3]/4)+(pi*radii1[3]*radii1[4]/4)+(pi*radii1[1]*radii1[4]/4)#
		#print(ba1)
		ba2<-(pi*radii2[1]*radii2[2]/4)+(pi*radii2[2]*radii2[3]/4)+(pi*radii2[3]*radii2[4]/4)+(pi*radii2[1]*radii2[4]/4)
		#print(ba2)
		bai<-ba1-ba2
		rgr<-(ba1-ba2)/ba1#relative growth rate=bai/ba
		#print(bai)
		indbai<-c(indbai,bai)#tacks on next) oldest year of bai to end of indbai
		indrgr<-c(indrgr,rgr)#tacks on next) oldest year of rgr to end of indrgr
		}
		#bai data
	if(lastring==2007){indbai2<-c(NA,NA,NA,indbai)}
	if(lastring==2008){indbai2<-c(NA,NA,indbai)}
	if(lastring==2009){indbai2<-c(NA,indbai)}
	if(lastring==2010){indbai2<-indbai}
		indbai3<-c(indbai2,rep(NA,times=(101-(length(indbai2)))))
		baidat<-rbind(baidat,indbai3)
		indsapdat<-tmpdat[1,5:16]
		inddat<-rbind(inddat,indsapdat)
		#rgr data
		if(lastring==2007){indrgr2<-c(NA,NA,NA,indrgr)}
	if(lastring==2008){indrgr2<-c(NA,NA,indrgr)}
	if(lastring==2009){indrgr2<-c(NA,indrgr)}
	if(lastring==2010){indrgr2<-indrgr}
		indrgr3<-c(indrgr2,rep(NA,times=(101-(length(indrgr2)))))
		rgrdat<-rbind(rgrdat,indrgr3)
	}
		colnames(baidat)<-seq(from=2010, to=1910)
		colnames(rgrdat)<-seq(from=2010, to=1910)
		sapbaidat<-cbind(inddat,baidat)
		saprgrdat<-cbind(inddat,rgrdat)
		
sapbaidat[1:5,1:17]
saprgrdat[1:5,1:17]
#write.csv(sapbaidat,"sapbaidat.csv")
#write.csv(saprgrdat,"saprgrdat.csv")
#If starting here, with baidat csv:
#sapbaidat<-read.csv("sapbaidat.csv",header=T)
#sapbaidat<-sapbaidat[,-1]#remove weird first column

###Climate data:
#Kevin's climate data estimates for each stand, which are based on PRISM estimates
standtemp<-read.csv("tree_plot_climate_temp.csv", header=TRUE)
standprecip<-read.csv("tree_plot_climate_precip.csv", header=TRUE)
unique(standprecip[,"Plot"])
dimnames(standprecip)
stands<-standtemp[,1]
month<-rep(1:12,length.out=1204)
year<-rep(1909:2009, each=12,length.out=1204)
monthyear<-cbind(year,month)
sttemp<-standtemp[,5:1208]
sttemp<-t(sttemp)
stdtemp<-cbind(monthyear,sttemp)
stdtemp<-stdtemp[61:1186,]#because snow data only goes back to 1914 and tree ring data ends at 2007 for most of the trees
colnames(stdtemp)<-c("year","month","AB08","AE10","AG05","AM16","AO03","AR07","AV02","AV06","AV14","AX15","PARA","PP17","TA01","TB13","TO04","TO11")
stprecip<-standprecip[,5:1208]
stprecip<-t(stprecip)
stdprecip<-cbind(monthyear,stprecip)
stdprecip<-stdprecip[61:1186,]
colnames(stdprecip)<-c("year","month","AB08","AE10","AG05","AM16","AO03","AR07","AV02","AV06","AV14","AX15","PARA","PP17","TA01","TB13","TO04","TO11")
row.names(stdtemp)<-NULL
row.names(stdprecip)<-NULL
sites<-c("PARA","AE10","AR07","AM16","AX15","AV06","AG05","TB13","TO04")
stdswe<-read.csv("tree_plot_climate_swe.csv", header=TRUE)
stdsnowdur<-read.csv("tree_plot_climate_snowdur.csv", header=TRUE)
stdgrdd5<-read.csv("tree_plot_climate_grdd5.csv", header=TRUE)
dim(stdsnowdur)
dim(stdswe)
dim(stdtemp)

#create a for loop to get all climate variables for each stand and examine correlations between them
for (i in 1:length(sites)){
	sitetemp<-cbind(stdtemp[,1:2],stdtemp[,dimnames(stdtemp)[[2]]==sites[i]])
	mntemp<-tapply(sitetemp[,3],sitetemp[,1],mean)
	sitetempgr<-sitetemp[sitetemp[,2]>5&sitetemp[,2]<11,]
	growtemp<-tapply(sitetempgr[,3],sitetempgr[,1],mean)
	sitetempdr<-sitetemp[sitetemp[,2]<5,]
	sitetempdr2<-sitetemp[sitetemp[,2]==12,]
	sitetempdr2[,1]<-sitetempdr2[,1]+1
	sitetempdr<-rbind(sitetempdr,sitetempdr2)
	dortemp<-tapply(sitetempdr[,3],sitetempdr[,1],mean)
	#now precip
	siteprecip<-cbind(stdprecip[,1:2],stdprecip[,dimnames(stdprecip)[[2]]==sites[i]])
	totprecip<-tapply(siteprecip[,3],siteprecip[,1],sum)
	siteprecipgr<-siteprecip[siteprecip[,2]>5&siteprecip[,2]<11,]
	growprecip<-tapply(siteprecipgr[,3],siteprecipgr[,1],sum)
	siteprecipdr<-siteprecip[siteprecip[,2]<5,]
	siteprecipdr2<-siteprecip[siteprecip[,2]==12,]
	siteprecipdr2[,1]<-siteprecipdr2[,1]+1
	siteprecipdr<-rbind(siteprecipdr,siteprecipdr2)
	dorprecip<-tapply(siteprecipdr[,3],siteprecipdr[,1],sum)
	#now snowdata&growingdegreedays(5c threshold)
	siteswe<-stdswe[,dimnames(stdswe)[[2]]==sites[i]]
	names(siteswe)<-stdswe[,1]
	swe<-siteswe[1:94]
	sitesnowdur<-stdsnowdur[,dimnames(stdsnowdur)[[2]]==sites[i]]
	names(sitesnowdur)<-stdsnowdur[,1]
	snowdur<-sitesnowdur[1:94]
	sitegrowdd<-stdgrdd5[,dimnames(stdgrdd5)[[2]]==sites[i]]
	names(sitegrowdd)<-stdgrdd5[,1]
	growdd<-sitegrowdd[1:94]
	#create matrix with all climate variables clim using cbind
	climsite<-cbind(mntemp,growtemp,dortemp,totprecip,growprecip,dorprecip,swe,snowdur,growdd)
	climsite<-climsite[order(dimnames(climsite)[[1]],decreasing=TRUE),]#to reorder climsite so that it is in the same order as the tree ring data!
	#create matrices for each stand's climate data(unstandardized)
	if(sites[i]=="PARA"){climPARA<-climsite}
	if(sites[i]=="AE10"){climAE10<-climsite}
	if(sites[i]=="AR07"){climAR07<-climsite}
	if(sites[i]=="AV06"){climAV06<-climsite}
	if(sites[i]=="AM16"){climAM16<-climsite}
	if(sites[i]=="AX15"){climAX15<-climsite}
	if(sites[i]=="AG05"){climAG05<-climsite}
	if(sites[i]=="TB13"){climTB13<-climsite}
	if(sites[i]=="TO04"){climTO04<-climsite}
	}

#Fit models of tree and sapling growth~climate  and growth~competition at each elevation. Then extract coefficients from models and plot coefficients vs elevation.
#sapbaidat: contains sapling bai data
#saprgrdat: contains sapling rgr data
#treebaidat: contains adult tree bai data
#treergrdat: contains adult tree rgr data

##First BAI
sapspecies<-unique(sapbaidat$Species)
stands<-c("TO04", "TB13","AG05","AV06", "AX15","AM16","AR07","AE10","PARA")
standelevs<-c(704,851,950,1064,1091,1197,1454,1460,1603)
allspmodout<-c()	
for (i in 1: length(species)){
	spdat<-sapbaidat[sapbaidat$Species==species[i],]
	allmodout<-c()
	for(j in 1:length(stands)){
		
 		spstanddat<-spdat[spdat$Stand.1==stands[j],]
 			if(dim(spstanddat)[1]==0){next}
			if(stands[j]=="PARA"){clim<-climPARA}
			if(stands[j]=="AE10"){clim<-climAE10}
			if(stands[j]=="AR07"){clim<-climAR07}
			if(stands[j]=="AV06"){clim<-climAV06}
			if(stands[j]=="AM16"){clim<-climAM16}
			if(stands[j]=="AX15"){clim<-climAX15}
			if(stands[j]=="AG05"){clim<-climAG05}
			if(stands[j]=="TB13"){clim<-climTB13}
			if(stands[j]=="TO04"){clim<-climTO04}
			recentMAT<-clim[1:12,2]#1996-2007 climate (MAT)
			recentSNDR<-clim[1:12,8]
			standdat<-c() 
			ind<-unique(spstanddat$TreeName)
			for(k in 1:length(ind)){
				inddat<-as.numeric(spstanddat[spstanddat	$TreeName==ind[k],16:27])#growth 1996-2007
 				#if(dim(inddat)[1]==0){next}
 				yrs<-seq(2007,1996)
			inds<-rep(k, times=length(yrs))
			gap<-rep(spstanddat[spstanddat				$TreeName==ind[k],6], times=length(yrs))
			tmp<-cbind(yrs,inds,gap,inddat)
				#put data into standdat (for later analysis)
			tmp3<-cbind(yrs,recentMAT,recentSNDR,tmp)
			#tmp4<-tmp3[is.na(tmp3[,11])==F,]
			standdat<-rbind(standdat,tmp3)
			}
			standdat[,1]<-as.factor(standdat[,1])
			standdat[,2]<-as.numeric(standdat[,2])
			standdat[,3]<-as.numeric(standdat[,3])
			
			standdat[,5]<-as.factor(standdat[,5])
			standdat[,6]<-as.factor(standdat[,6])
			standdat[,7]<-as.numeric(standdat[,7])
			#models
			MAT<-as.numeric(standdat[,2])
			SNDR<-as.numeric(standdat[,3])
			year<-as.factor(standdat[,4])
			inds<-as.factor(standdat[,5])
			gap<-as.factor(standdat[,6])
			growth<-as.numeric(standdat[,7])
			matmod<-lm(log(growth)~MAT)#growth needs to be logged to be normal
			snwmod<-lm(log(growth)~SNDR)#growth needs to be logged to be normal
			compmod<-lm(log(growth)~gap)
			compcoef<-summary(compmod)[[4]][2,1]#gap effect
			compcoefp<-summary(compmod)[[4]][2,4]#pval gap
			compR2<-summary(compmod)[9]#adj R2 comp. mod
			compp<-summary(compmod)[[4]][2,4]#p-val comp.mod
			snwcoef<-summary(snwmod)[[4]][2,1]#SNWDR effect
			snwcoefp<-summary(snwmod)[[4]][2,4]#pval SNDR
			snwR2<-summary(snwmod)[9]#adj R2 snwdur mod
			snwp<-summary(snwmod)[[4]][2,4]#p-val snwmod
			spec<-as.character(sapspecies[i])
			elev<-standelevs[j]
			modout<-cbind(spec,elev,compcoef,compcoefp,compR2,compp,snwcoef,snwcoefp,snwR2,snwp)
			allmodout<-rbind(allmodout,modout)
			}
			allspmodout<-rbind(allspmodout,allmodout)
			}
colnames(allspmodout)<-c("Species","Elev","CompModCoef","CompModCoefP","CompModR2","CompModP","SNDRModCoef","SNDRModCoefP","SNDRModR2","SNDRModP")
write.csv(allspmodout,"sapmodelsbai.csv")

#Now adult trees
species<-c("Abam","Tshe","Tsme")
stands<-c("TO04", "TB13","AG05","AV06", "AX15","AM16","AR07","AE10","PARA")
standelevs<-c(704,851,950,1064,1091,1197,1454,1460,1603)
alladspmodout<-c()	
for (i in 1: length(species)){
	adspdat<-treebaidat[treebaidat[,3]==species[i],]
	alladmodout<-c()
	for(j in 1:length(stands)){
 		adspstanddat<-adspdat[adspdat[,1]==stands[j],]
 			if(dim(adspstanddat)[1]==0){next}
			if(stands[j]=="PARA"){clim<-climPARA}
			if(stands[j]=="AE10"){clim<-climAE10}
			if(stands[j]=="AR07"){clim<-climAR07}
			if(stands[j]=="AV06"){clim<-climAV06}
			if(stands[j]=="AM16"){clim<-climAM16}
			if(stands[j]=="AX15"){clim<-climAX15}
			if(stands[j]=="AG05"){clim<-climAG05}
			if(stands[j]=="TB13"){clim<-climTB13}
			if(stands[j]=="TO04"){clim<-climTO04}
			recentMAT<-clim[1:12,2]#1996-2007 climate (MAT)
			recentSNDR<-clim[1:12,8]
			standdat<-c() 
			ind<-unique(adspstanddat[,4])
			for(k in 1:length(ind)){
				inddat<-as.numeric(adspstanddat[adspstanddat	[,4]==ind[k],7:18])#growth 1996-2007
 				#if(dim(inddat)[1]==0){next}
 				yrs<-seq(2007,1996)
			inds<-rep(k, times=length(yrs))
			#gap<-rep(spstanddat[spstanddat				[,4]==ind[k],6], times=length(yrs))#need to add gap data
			tmp<-cbind(yrs,inds,inddat)#add "gap" after inds, once you calculate it
				#put data into standdat (for later analysis)
			tmp3<-cbind(yrs,recentMAT,recentSNDR,tmp)
			#tmp4<-tmp3[is.na(tmp3[,11])==F,]
			standdat<-rbind(standdat,tmp3)
			}
			standdat[,1]<-as.factor(standdat[,1])
			standdat[,2]<-as.numeric(standdat[,2])
			standdat[,3]<-as.numeric(standdat[,3])
			standdat[,4]<-as.factor(standdat[,4])
			standdat[,5]<-as.factor(standdat[,5])
			standdat[,6]<-as.numeric(standdat[,6])
			#models
			MAT<-as.numeric(standdat[,2])
			SNDR<-as.numeric(standdat[,3])
			year<-as.factor(standdat[,4])
			inds<-as.factor(standdat[,5])
			#gap<-as.factor(standdat[,6])
			growth<-as.numeric(standdat[,6])
			matmod<-lm(log(growth)~MAT)#growth needs to be logged to be normal
			print(summary(matmod))
			snwmod<-lm(log(growth)~SNDR)#growth needs to be logged to be normal
			print(summary(snwmod))
			#compmod<-lm(log(growth)~gap)
			#compcoef<-summary(compmod)[[4]][2,1]#gap effect
			#compcoefp<-summary(compmod)[[4]][2,4]#pval gap
			#compR2<-summary(compmod)[9]#adj R2 comp. mod
			#compp<-summary(compmod)[[4]][2,4]#p-val comp.mod
			snwcoef<-summary(snwmod)[[4]][2,1]#SNWDR effect
			snwcoefp<-summary(snwmod)[[4]][2,4]#pval SNDR
			snwR2<-summary(snwmod)[9]#adj R2 snwdur mod
			snwp<-summary(snwmod)[[4]][2,4]#p-val snwmod
			spec<-as.character(sapspecies[i])
			elev<-standelevs[j]
			#modout<-cbind(spec,elev,compcoef,compcoefp,compR2,compp,snwcoef,snwcoefp,snwR2,snwp)
			modout<-cbind(spec,elev,snwcoef,snwcoefp,snwR2,snwp)#add comp stuff later
			alladmodout<-rbind(alladmodout,modout)
			}
			alladspmodout<-rbind(alladspmodout,alladmodout)
			}
#colnames(alladspmodout)<-c("Species","Elev","CompModCoef","CompModCoefP","CompModR2","CompModP","SNDRModCoef","SNDRModCoefP","SNDRModR2","SNDRModP")
colnames(alladspmodout)<-c("Species","Elev","SNDRModCoef","SNDRModCoefP","SNDRModR2","SNDRModP")#add comp stuff later
write.csv(alladspmodout,"adultmodelsbai.csv")
#Plot for ABAM BAI
quartz(height=6, width=10)
par(mfrow=c(2,2))
#saplings & competition
compfill<-rep("white", length=7)
compfill[allspmodout[1:7,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[1:7,3]), names.arg=allspmodout[1:7,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#Adults & competition
barplot(rep(0,times=8), names.arg=alladspmodout[1:8,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=7)
snwfill[allspmodout[1:7,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[1:7,7]), names.arg=allspmodout[1:7,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill, ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=8)
adsnwfill[alladspmodout[1:8,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[1:8,3]), names.arg=alladspmodout[1:8,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))

#Plot for TSHE BAI
quartz(height=6, width=10)
par(mfrow=c(2,2))
#saplings & competition
compfill<-rep("white", length=4)
compfill[allspmodout[8:11,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[8:11,3]), names.arg=allspmodout[8:11,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#adults & competition
barplot(rep(0,times=6), names.arg=alladspmodout[9:14,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=4)
snwfill[allspmodout[8:11,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[8:11,7]), names.arg=allspmodout[8:11,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill,ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=6)
adsnwfill[alladspmodout[9:14,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[9:14,3]), names.arg=alladspmodout[9:14,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))


#Plot for TSME BAI
quartz(height=6, width=10)
par(mfrow=c(2,2))
compfill<-rep("white", length=4)
compfill[allspmodout[12:15,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[12:15,3]), names.arg=allspmodout[12:15,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#adults & competition
barplot(rep(0,times=4), names.arg=alladspmodout[15:18,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=4)
snwfill[allspmodout[12:15,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[12:15,7]), names.arg=allspmodout[12:15,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill, ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=4)
adsnwfill[alladspmodout[15:18,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[15:18,3]), names.arg=alladspmodout[15:18,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))



##Next RGR
#saprgrdat: contains sapling rgr data
#treergrdat: contains adult tree rgr data
allspmodout<-c()	
for (i in 1: length(species)){
	spdat<-saprgrdat[saprgrdat$Species==species[i],]
	allmodout<-c()
	for(j in 1:length(stands)){
		
 		spstanddat<-spdat[spdat$Stand.1==stands[j],]
 			if(dim(spstanddat)[1]==0){next}
			if(stands[j]=="PARA"){clim<-climPARA}
			if(stands[j]=="AE10"){clim<-climAE10}
			if(stands[j]=="AR07"){clim<-climAR07}
			if(stands[j]=="AV06"){clim<-climAV06}
			if(stands[j]=="AM16"){clim<-climAM16}
			if(stands[j]=="AX15"){clim<-climAX15}
			if(stands[j]=="AG05"){clim<-climAG05}
			if(stands[j]=="TB13"){clim<-climTB13}
			if(stands[j]=="TO04"){clim<-climTO04}
			recentMAT<-clim[1:12,2]#1996-2007 climate (MAT)
			recentSNDR<-clim[1:12,8]
			standdat<-c() 
			ind<-unique(spstanddat$TreeName)
			for(k in 1:length(ind)){
				inddat<-as.numeric(spstanddat[spstanddat	$TreeName==ind[k],16:27])#growth 1996-2007
 				#if(dim(inddat)[1]==0){next}
 				yrs<-seq(2007,1996)
			inds<-rep(k, times=length(yrs))
			gap<-rep(spstanddat[spstanddat				$TreeName==ind[k],6], times=length(yrs))
			tmp<-cbind(yrs,inds,gap,inddat)
				#put data into standdat (for later analysis)
			tmp3<-cbind(yrs,recentMAT,recentSNDR,tmp)
			#tmp4<-tmp3[is.na(tmp3[,11])==F,]
			standdat<-rbind(standdat,tmp3)
			}
			standdat[,1]<-as.factor(standdat[,1])
			standdat[,2]<-as.numeric(standdat[,2])
			standdat[,3]<-as.numeric(standdat[,3])
			
			standdat[,5]<-as.factor(standdat[,5])
			standdat[,6]<-as.factor(standdat[,6])
			standdat[,7]<-as.numeric(standdat[,7])
			#models
			MAT<-as.numeric(standdat[,2])
			SNDR<-as.numeric(standdat[,3])
			year<-as.factor(standdat[,4])
			inds<-as.factor(standdat[,5])
			gap<-as.factor(standdat[,6])
			growth<-as.numeric(standdat[,7])
			matmod<-lm(log(growth)~MAT)#growth needs to be logged to be normal
			snwmod<-lm(log(growth)~SNDR)#growth needs to be logged to be normal
			compmod<-lm(log(growth)~gap)
			compcoef<-summary(compmod)[[4]][2,1]#gap effect
			compcoefp<-summary(compmod)[[4]][2,4]#pval gap
			compR2<-summary(compmod)[9]#adj R2 comp. mod
			compp<-summary(compmod)[[4]][2,4]#p-val comp.mod
			snwcoef<-summary(snwmod)[[4]][2,1]#SNWDR effect
			snwcoefp<-summary(snwmod)[[4]][2,4]#pval SNDR
			snwR2<-summary(snwmod)[9]#adj R2 snwdur mod
			snwp<-summary(snwmod)[[4]][2,4]#p-val snwmod
			spec<-as.character(sapspecies[i])
			elev<-standelevs[j]
			modout<-cbind(spec,elev,compcoef,compcoefp,compR2,compp,snwcoef,snwcoefp,snwR2,snwp)
			allmodout<-rbind(allmodout,modout)
			}
			allspmodout<-rbind(allspmodout,allmodout)
			}
colnames(allspmodout)<-c("Species","Elev","CompModCoef","CompModCoefP","CompModR2","CompModP","SNDRModCoef","SNDRModCoefP","SNDRModR2","SNDRModP")


#Now adult trees
species<-c("Abam","Tshe","Tsme")
stands<-c("TO04", "TB13","AG05","AV06", "AX15","AM16","AR07","AE10","PARA")
standelevs<-c(704,851,950,1064,1091,1197,1454,1460,1603)
alladspmodout<-c()	
for (i in 1: length(species)){
	adspdat<-treergrdat[treergrdat[,3]==species[i],]
	alladmodout<-c()
	for(j in 1:length(stands)){
 		adspstanddat<-adspdat[adspdat[,1]==stands[j],]
 			if(dim(adspstanddat)[1]==0){next}
			if(stands[j]=="PARA"){clim<-climPARA}
			if(stands[j]=="AE10"){clim<-climAE10}
			if(stands[j]=="AR07"){clim<-climAR07}
			if(stands[j]=="AV06"){clim<-climAV06}
			if(stands[j]=="AM16"){clim<-climAM16}
			if(stands[j]=="AX15"){clim<-climAX15}
			if(stands[j]=="AG05"){clim<-climAG05}
			if(stands[j]=="TB13"){clim<-climTB13}
			if(stands[j]=="TO04"){clim<-climTO04}
			recentMAT<-clim[1:12,2]#1996-2007 climate (MAT)
			recentSWE<-clim[1:12,7]
			standdat<-c() 
			ind<-unique(adspstanddat[,4])
			for(k in 1:length(ind)){
				inddat<-as.numeric(adspstanddat[adspstanddat	[,4]==ind[k],7:18])#growth 1996-2007
 				#if(dim(inddat)[1]==0){next}
 				yrs<-seq(2007,1996)
			inds<-rep(k, times=length(yrs))
			#gap<-rep(spstanddat[spstanddat				[,4]==ind[k],6], times=length(yrs))#need to add gap data
			tmp<-cbind(yrs,inds,inddat)#add "gap" after inds, once you calculate it
				#put data into standdat (for later analysis)
			tmp3<-cbind(yrs,recentMAT,recentSWE,tmp)
			#tmp4<-tmp3[is.na(tmp3[,11])==F,]
			standdat<-rbind(standdat,tmp3)
			}
			standdat[,1]<-as.factor(standdat[,1])
			standdat[,2]<-as.numeric(standdat[,2])
			standdat[,3]<-as.numeric(standdat[,3])
			standdat[,4]<-as.factor(standdat[,4])
			standdat[,5]<-as.factor(standdat[,5])
			standdat[,6]<-as.numeric(standdat[,6])
			#models
			MAT<-as.numeric(standdat[,2])
			SWE<-as.numeric(standdat[,3])
			year<-as.factor(standdat[,4])
			inds<-as.factor(standdat[,5])
			#gap<-as.factor(standdat[,6])
			growth<-as.numeric(standdat[,6])
			matmod<-lm(log(growth)~MAT)#
			print(summary(matmod))
			snwmod<-lm(log(growth)~SWE)#growth needs to be logged to be normal
			print(summary(snwmod))
			#compmod<-lm(log(growth)~gap)
			#compcoef<-summary(compmod)[[4]][2,1]#gap effect
			#compcoefp<-summary(compmod)[[4]][2,4]#pval gap
			#compR2<-summary(compmod)[9]#adj R2 comp. mod
			#compp<-summary(compmod)[[4]][2,4]#p-val comp.mod
			snwcoef<-summary(snwmod)[[4]][2,1]#SWE effect
			snwcoefp<-summary(snwmod)[[4]][2,4]#pval SWE
			snwR2<-summary(snwmod)[9]#adj R2 snw mod
			snwp<-summary(snwmod)[[4]][2,4]#p-val snwmod
			spec<-as.character(sapspecies[i])
			elev<-standelevs[j]
			#modout<-cbind(spec,elev,compcoef,compcoefp,compR2,compp,snwcoef,snwcoefp,snwR2,snwp)
			modout<-cbind(spec,elev,snwcoef,snwcoefp,snwR2,snwp)#add comp stuff later
			alladmodout<-rbind(alladmodout,modout)
			}
			alladspmodout<-rbind(alladspmodout,alladmodout)
			}
#colnames(alladspmodout)<-c("Species","Elev","CompModCoef","CompModCoefP","CompModR2","CompModP","SNDRModCoef","SNDRModCoefP","SNDRModR2","SNDRModP")
colnames(alladspmodout)<-c("Species","Elev","SWEModCoef","SWEModCoefP","SWEModR2","SWEModP")#add comp stuff later





#Plot for ABAM RGR
quartz(height=6, width=10)
par(mfrow=c(2,2))
#saplings & competition
compfill<-rep("white", length=7)
compfill[allspmodout[1:7,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[1:7,3]), names.arg=allspmodout[1:7,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#Adults & competition
barplot(rep(0,times=8), names.arg=alladspmodout[1:8,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=7)
snwfill[allspmodout[1:7,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[1:7,7]), names.arg=allspmodout[1:7,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill, ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=8)
adsnwfill[alladspmodout[1:8,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[1:8,3]), names.arg=alladspmodout[1:8,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))

#Plot for TSHE RGR
quartz(height=6, width=10)
par(mfrow=c(2,2))
#saplings & competition
compfill<-rep("white", length=4)
compfill[allspmodout[8:11,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[8:11,3]), names.arg=allspmodout[8:11,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#adults & competition
barplot(rep(0,times=6), names.arg=alladspmodout[9:14,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=4)
snwfill[allspmodout[8:11,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[8:11,7]), names.arg=allspmodout[8:11,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill, ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=6)
adsnwfill[alladspmodout[9:14,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[9:14,3]), names.arg=alladspmodout[9:14,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))


#Plot for TSME RGR
quartz(height=6, width=10)
par(mfrow=c(2,2))
compfill<-rep("white", length=4)
compfill[allspmodout[12:15,4]<0.05]<-"black"
barplot(as.numeric(allspmodout[12:15,3]), names.arg=allspmodout[12:15,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of competition on growth", col=compfill, main="Sapling Growth")

#adults & competition
barplot(rep(0,times=4), names.arg=alladspmodout[15:18,2], beside=TRUE, xlab="Elevation (m)", col=compfill, main="Adult Tree Growth")

#saplings & climate
snwfill<-rep("white", length=4)
snwfill[allspmodout[12:15,8]<0.05]<-"black"
barplot(as.numeric(allspmodout[12:15,7]), names.arg=allspmodout[12:15,2], beside=TRUE, xlab="Elevation (m)", ylab="Effect of snow duration on growth", col=snwfill,ylim=c(-0.005,0.005))

#adults & climate
adsnwfill<-rep("white", length=4)
adsnwfill[alladspmodout[15:18,4]<0.05]<-"black"
barplot(as.numeric(alladspmodout[15:18,3]), names.arg=alladspmodout[15:18,2], beside=TRUE, xlab="Elevation (m)", col=adsnwfill,ylim=c(-0.005,0.005))


###use basal area of neighbors to create categorical variable for competition (like gap/nongap) for adult trees
#first, read in tree competive/neighbor data
treenbrs<-read.csv("SouthSideNbrsUpdate.csv", header=TRUE)
head(treenbrs)
dim(treenbrs)
unique(treenbrs[,3])#460 focal trees/plots
elev<-as.numeric(treenbrs$Elevation)#elevation
dist<-as.numeric(treenbrs$dist)#Distance
dbh<-as.numeric(treenbrs$Dbh)#Dbh
treeid<-as.factor(treenbrs$Tree)#Tag number
stands<-unique(treenbrs$Stand)
basarea<-pi*(dbh/2)*(dbh/2)
banbrs<-tapply(basarea,treeid,sum, na.rm = TRUE)#total basal area of neighbors around each focal tree- a measure of competition
#sapling competitive/neighbor data
sapnbrs<-read.csv("SaplingNbrsSouthSide.csv",header=T)
unique(sapnbrs$New.)#316 unique saplings in nbr data
#first, get basal area of each neighbot
basarea<-pi*(as.numeric(sapnbrs$Neighbor.DBH)/2)*(as.numeric(sapnbrs$Neighbor.DBH)/2)
comp<-tapply(basarea,sapnbrs$New.,sum, na.rm = TRUE)
max(comp)
which(names(comp=="202827.0"))
#area basal areas of adult trees and saplings tree neighbors similar?
par(mfrow=c(1,2))
hist(banbrs)
hist(comp[2:316])#remove first entry, which doesn't go with a sapling...not sure why???
saplings<-rep("saplings", times=length(comp[2:316]))
adults<-rep("adults", times=length(banbrs))
stage<-c(adults,saplings)
ba<-c(banbrs,comp[2:316])
t.test(ba~stage)

#look similar based on histogram, but t.test says that means are different:Welch Two Sample t-test data:  ba by stage t = 1.9624, df = 662.921, p-value =0.05013, alternative hypothesis: true difference in means is not equal to 0 , 95 percent confidence interval:-0.9668817 3409.8102972, sample estimates:, mean in group adults=27069.70,  mean in group saplings= 25365.27 


#Look at total basal area of neighbors in gaps vs. nongaps at each stand, 5 m from focal tree. Use this to create a categorical variable for competitive environment for adult trees. 
dist5<-sapnbrs[sapnbrs$Distance<5,]#just count neighbors within 5 m
badist5<-pi*(as.numeric(dist5$Neighbor.DBH)/2)*(as.numeric(dist5$Neighbor.DBH)/2)
allbadist5<-tapply(badist5,dist5$New.,sum,na.rm = TRUE)
allbadist5[which(is.na(allbadist5))]<-0
gapstat<-c()
for(i in 1:length(sapbaidat[,1])){
		indgapstat<-sapbaidat[i,6]#gap status of tree (y/n)
		sapcomp<-allbadist5[names(allbadist5)==sapbaidat$TreeName[i]]#basal area for neighbors of focal tree
		std<-sapbaidat$Stand.1[i]
		sp<-sapbaidat$Species[i]
		if(length(sapcomp)[1]==0){next}
		indgapstatcomp<-c(as.character(sp),as.character(std),names(sapcomp),as.character(indgapstat),as.numeric(sapcomp))
		gapstat<-rbind(gapstat,indgapstatcomp)
		}
		colnames(gapstat)<-c("Species","Stand","TreeName","gapstat","compba5")
		#gapstatall<-as.factor(gapstat[,1])

head(gapstat)

basarea<-as.numeric(gapstat[,5])
gap<-gapstat[,4]
gap[which(gap=="y")]<-1
gap[which(gap=="n")]<-0
gap<-as.numeric(gap)
gapmod<-glm(gap~basarea,family=binomial)
summary(gapmod)
plot(gap~basarea)
basarea<-seq(min(basarea),max(basarea),length.out=273)
gap.predicted<-predict.glm(gapmod,as.data.frame(basarea),type="response")

lines(basarea,gap.predicted)
basarea<-as.numeric(gapstat[,5])

###Exploring the data a bit, for saplings:
#Plot RGR~stand and BAI~stand,for all 3 species to see if there are stand differences in growth of saplings:
abamstands<-c("PARA","AE10","AR07","AM16","AV06","TB13","TO04")#ordered from high to low elevation
#unique(tshergrdat$Stand.1)
tsmestands<-c("PARA","AE10","AR07","AM16")
tshestands<-c("AM16","AV06","TB13","TO04")
abamrgrdat<-saprgrdat[saprgrdat$Species=="Abam",]
#Now try averaging growth across most recent 12 years	
abammngrwth<-rowMeans(abamrgrdat[,15:27], na.rm = TRUE)
tsmergrdat<-saprgrdat[saprgrdat$Species=="Tsme",]
tsmemngrwth<-rowMeans(tsmergrdat[,15:27], na.rm = TRUE)
tshergrdat<-saprgrdat[saprgrdat$Species=="Tshe",]
tshemngrwth<-rowMeans(tshergrdat[,15:27], na.rm = TRUE)
abambaidat<-sapbaidat[sapbaidat$Species=="Abam",]
abammnbai<-rowMeans(abambaidat[,15:27], na.rm = TRUE)
tsmebaidat<-sapbaidat[sapbaidat$Species=="Tsme",]
tsmemnbai<-rowMeans(tsmebaidat[,15:27], na.rm = TRUE)
tshebaidat<-sapbaidat[sapbaidat$Species=="Tshe",]
tshemnbai<-rowMeans(tshebaidat[,15:27], na.rm = TRUE)
#RGR
quartz(width=12,height=5)
par(mfrow=c(1,3))
boxplot(abammngrwth~abamrgrdat$Stand.1, xlab="Stand", ylab="Recent Mean RGR", main="ABAM")
boxplot(tshemngrwth~tshergrdat$Stand.1, xlab="Stand", ylab="Recent Mean RGR", main="TSHE")
boxplot(tsmemngrwth~tsmergrdat$Stand.1, xlab="Stand", ylab="Recent Mean RGR",main="TSME")
summary(aov(abammngrwth~abamrgrdat$Stand.1))#differences are significant (P<0.001)
summary(aov(tshemngrwth~tshergrdat$Stand.1))#differences are significant (P<0.001)
summary(aov(tsmemngrwth~tsmergrdat$Stand.1))#no significant differences (P=0.58)
#BAI
quartz(width=12,height=5)
par(mfrow=c(1,3))
boxplot(abammnbai~abambaidat$Stand.1, xlab="Stand", ylab="Recent Mean BAI", main="ABAM")
boxplot(tshemnbai~tshebaidat$Stand.1, xlab="Stand", ylab="Recent Mean BAI", main="TSHE")
boxplot(tsmemnbai~tsmebaidat$Stand.1, xlab="Stand", ylab="Recent Mean BAI",main="TSME")
summary(aov(abammnbai~abambaidat$Stand.1))#differences are not significant (P=0.38)
summary(aov(tshemnbai~tshebaidat$Stand.1))#differences are marginal (P=0.096)
summary(aov(tsmemnbai~tsmebaidat$Stand.1))#no significant differences (P=0.214)

#RGR and BAI in gaps vs nongaps
#ABAM
#rgr
abamrgrdat<-saprgrdat[saprgrdat$Species=="Abam",]
quartz(height=7,width=11)
par(mfrow=c(3,3))
for(i in 1:length(abamstands)){
	tmpdat<-abamrgrdat[abamrgrdat$Stand.1==abamstands[i],]
	stdgrthopen<-c()
	for (j in 1:dim(tmpdat)[1]){
		mnsapgr<-rowMeans(tmpdat[j,16:27],na.rm = TRUE)#mean bai frmo 1996-2007 years
		popen<-tmpdat[j,5]#densiometer reading
		grthopen<-cbind(mnsapgr,popen)
		stdgrthopen<-rbind(stdgrthopen,grthopen)
		}
		colnames(stdgrthopen)<-c("mngrwth","pcntopen")
		rownames(stdgrthopen)<-tmpdat[,2]
		openpcnt<-stdgrthopen[,2]
		mngrth<-stdgrthopen[,1]
plot(openpcnt,mngrth,xlab="Percent Open",ylab="",ylim=c(0,0.3), main=paste("Abam Sapling Growth",abamstands[i]), pch=16,cex=1.5, cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
if (abamstands[i]=="AM16"){mtext("Average annual growth (BAI, mm2)", side=2, line=2,cex=1.2)}
mod<-lm(mngrth~openpcnt)
modresults<-summary(mod)
print(modresults)
abline(mod,lwd=3)
mtext(paste("p=",round(modresults$coef[8], digits=3)), side=1, at=.9*max(openpcnt), cex=.9, line=-0.8)
mtext(paste("R2=",round(modresults$adj.r.squared,digits=3)),side=1, at=.9*max(openpcnt), cex=.9, line=-1.8)
}

#bai
quartz(height=7,width=11)
par(mfrow=c(3,3))
for(i in 1:length(stands)){
	tmpdat<-abambaidat[abambaidat$Stand.1==abamstands[i],]
	stdgrthopen<-c()
	for (j in 1:dim(tmpdat)[1]){
		mnsapgr<-rowMeans(tmpdat[j,16:27],na.rm = TRUE)#mean bai frmo 1996-2007 years
		popen<-tmpdat[j,5]#densiometer reading
		grthopen<-cbind(mnsapgr,popen)
		stdgrthopen<-rbind(stdgrthopen,grthopen)
		}
		colnames(stdgrthopen)<-c("mngrwth","pcntopen")
		rownames(stdgrthopen)<-tmpdat[,2]
		openpcnt<-stdgrthopen[,2]
		mngrth<-stdgrthopen[,1]
plot(openpcnt,mngrth,xlab="Percent Open",ylab="",ylim=c(0,60), main=paste("Abam Sapling Growth",abamstands[i]), pch=16,cex=1.5, cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
if (abamstands[i]=="AM16"){mtext("Average annual growth (BAI, mm2)", side=2, line=2,cex=1.2)}
mod<-lm(mngrth~openpcnt)
modresults<-summary(mod)
print(modresults)
abline(mod,lwd=3)
mtext(paste("p=",round(modresults$coef[8], digits=3)), side=1, at=.9*max(openpcnt), cex=.9, line=-0.8)
mtext(paste("R2=",round(modresults$adj.r.squared,digits=3)),side=1, at=.9*max(openpcnt), cex=.9, line=-1.8)
}

#mean growth in gaps vs nongapsm using bai and RGR
abammnrgr<-rowMeans(abamrgrdat[,15:27], na.rm = TRUE)
abammeans<-tapply(abammnrgr,list(abamrgrdat$In.Gap.,abamrgrdat$Stand.1), mean, na.rm = TRUE)
abamsd<-tapply(abammnrgr,list(abamrgrdat$In.Gap.,abamrgrdat$Stand.1), sd, na.rm = TRUE)#
abamn<-tapply(abammnrgr,list(abamrgrdat$In.Gap.,abamrgrdat$Stand.1), length)
abamse<-abamsd/sqrt(abamn)
abamgaps<-abammeans[3,2:8]
abamnotgaps<-abammeans[2,2:8]
x<-c(6,4,5,3,7,2,1)
abamgaps<-abamgaps[order(x)]
abamnotgaps<-abamnotgaps[order(x)]
abamgapse<-abamse[3,2:8]
abamnotgapse<-abamse[2,2:8]
abamgapse<-abamgapse[order(x)]
abamnotgapse<-abamnotgapse[order(x)]
allabammeans<-c(abammeans[2:3,2:8])
allabamse<-c(abamse[2:3,2:8])
y<-c(12,11,8,7,10,9,6,5,14,13,4,3,2,1)
allabammeans<-allabammeans[order(y)]
allabamse<-allabamse[order(y)]

quartz(height=6,width=10)
par(mfrow=c(2,1))
abamplot<-barplot(as.matrix(rbind(abamgaps,abamnotgaps)),ylab="Mean RGR",xlab="Elevation (m)",ylim=c(0,0.15),las=1,col=c("darkorange3","dark blue"),beside=TRUE,names.arg=c("704","851","1064","1197","1454","1460","1600"), main="Abies amabilis sapling growth",cex.names=1.3,cex.lab=1.5,cex.main=1.8,cex.axis=1.3)
legend(13,0.15,legend=c("Gap","Non Gap"),bty="n",fill=c("darkorange3","dark blue"),cex=1.5)
x<-c(abamplot)
for (i in 1:length(allabammeans)){
	arrows(x[i],allabammeans[i]-allabamse[i],x[i],allabammeans[i]+allabamse[i],length=0.1,angle=90,code=3,lwd=2)}
#add BAI
abammnbai<-rowMeans(abambaidat[,15:27], na.rm = TRUE)
abammeans<-tapply(abammnbai,list(abambaidat$In.Gap.,abambaidat$Stand.1), mean, na.rm = TRUE)
abamsd<-tapply(abammnbai,list(abambaidat$In.Gap.,abambaidat$Stand.1), sd, na.rm = TRUE)#
abamn<-tapply(abammnbai,list(abambaidat$In.Gap.,abambaidat$Stand.1), length)
abamse<-abamsd/sqrt(abamn)
abamgaps<-abammeans[3,2:8]
abamnotgaps<-abammeans[2,2:8]
x<-c(6,4,5,3,7,2,1)
abamgaps<-abamgaps[order(x)]
abamnotgaps<-abamnotgaps[order(x)]
abamgapse<-abamse[3,2:8]
abamnotgapse<-abamse[2,2:8]
abamgapse<-abamgapse[order(x)]
abamnotgapse<-abamnotgapse[order(x)]
allabammeans<-c(abammeans[2:3,2:8])
allabamse<-c(abamse[2:3,2:8])
y<-c(12,11,8,7,10,9,6,5,14,13,4,3,2,1)
allabammeans<-allabammeans[order(y)]
allabamse<-allabamse[order(y)]
abamplot<-barplot(as.matrix(rbind(abamgaps,abamnotgaps)),ylab="Mean BAI",xlab="Elevation (m)",ylim=c(0,50),las=1,col=c("darkorange3","dark blue"),beside=TRUE,names.arg=c("704","851","1064","1197","1454","1460","1600"),cex.names=1.3,cex.lab=1.5,cex.axis=1.3)
x<-c(abamplot)
for (i in 1:length(allabammeans)){
	arrows(x[i],allabammeans[i]-allabamse[i],x[i],allabammeans[i]+allabamse[i],length=0.1,angle=90,code=3,lwd=2)}
