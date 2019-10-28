library(quantmod)
library(zoo)
# Gap, Inflation, FFR, GS2, GS5, GS10
for (series in c("UNRATE","NROU","TB3MS","GS2","GS5","GS10")) {
    getSymbols(series,src="FRED",.env=.GlobalEnv)
    eval(parse(text=paste0(series,"<-apply.quarterly(",series,",mean)")))
}
for (series in c("PCECTPI","GDPC1","GDPPOT","JCXFE")) {
    getSymbols(series,src="FRED",.env=.GlobalEnv)
}
ugap <- -(log(ts(UNRATE,start=1948,frequency=4)/100+1)-log(ts(NROU,start=1947,frequency=4)/100+1))*100
ygap <- log(ts(GDPC1,start=1947,frequency=4)/ts(GDPPOT,start=1949,frequency=4))*100
headline <- ts(log(PCECTPI/lag(PCECTPI,1))*400,start=1947,frequency=4)
core <- ts(log(JCXFE/lag(JCXFE,1))*400,start=1959,frequency=4)
r2 <- ts(log(GS2/100+1)*100,start=1976.25,frequency=4)
r5 <- ts(log(GS5/100+1)*100,start=1953.25,frequency=4)
r10 <- ts(log(GS10/100+1)*100,start=1953.25,frequency=4)
ffr <- ts(log(TB3MS/100+1)*100,start=1934,frequency=4)
ffr[ffr<0.25] <- NA

end.date <- 2018.75
start.date <- 1960

ugap.headline <- window(cbind(ugap,headline,ffr,r2,r5,r10),start=start.date,end=end.date)
ugap.core <- window(cbind(ugap,core,ffr,r2,r5,r10),start=start.date,end=end.date)
ygap.headline <- window(cbind(ygap,headline,ffr,r2,r5,r10),start=start.date,end=end.date)
ygap.core <- window(cbind(ygap,core,ffr,r2,r5,r10),start=start.date,end=end.date)

matlab.dates <- rep(NA,dim(ygap.headline)[1])
idx <- 0
for (y in 1960:2018) {
    for (m in c("01","04","07","10")) {
        idx <- idx + 1
	if (idx <= length(matlab.dates)) {
            matlab.dates[idx] <- 715876 + as.numeric(as.Date(paste0(y,"-",m,"-01"))) - as.numeric(as.Date("1960-01-01"))
        }
    }
}

prefix <- "spectreTB3MSGS020510"
datestr <- "2018Q4"

write.table(matrix(as.integer(is.na(ygap.headline)),dim(ygap.headline)[1]),file=paste0(prefix,"OutputGap","Headline",datestr,".yNaN.txt"),sep=" ",row.names=rep("",dim(ygap.headline)[1]),col.names=FALSE,quote=FALSE)
ygap.headline[is.na(ygap.headline)] <- 0
write.table(matrix(sprintf('%30.16E',ygap.headline),dim(ygap.headline)[1]),file=paste0(prefix,"OutputGap","Headline",datestr,".yData.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrix(sprintf('%30.16E',matlab.dates[1:dim(ygap.headline)[1]]),dim(ygap.headline)[1]),file=paste0(prefix,"OutputGap","Headline",datestr,".dates.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)


write.table(matrix(as.integer(is.na(ygap.core)),dim(ygap.core)[1]),file=paste0(prefix,"OutputGap","Core",datestr,".yNaN.txt"),sep=" ",row.names=rep("",dim(ygap.core)[1]),col.names=FALSE,quote=FALSE)
ygap.core[is.na(ygap.core)] <- 0
write.table(matrix(sprintf('%30.16E',ygap.core),dim(ygap.core)[1]),file=paste0(prefix,"OutputGap","Core",datestr,".yData.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrix(sprintf('%30.16E',matlab.dates[1:dim(ygap.core)[1]]),dim(ygap.core)[1]),file=paste0(prefix,"OutputGap","Core",datestr,".dates.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)


write.table(matrix(as.integer(is.na(ugap.core)),dim(ugap.core)[1]),file=paste0(prefix,"Core",datestr,".yNaN.txt"),sep=" ",row.names=rep("",dim(ugap.core)[1]),col.names=FALSE,quote=FALSE)
ugap.core[is.na(ugap.core)] <- 0
write.table(matrix(sprintf('%30.16E',ugap.core),dim(ugap.core)[1]),file=paste0(prefix,"Core",datestr,".yData.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrix(sprintf('%30.16E',matlab.dates[1:dim(ugap.core)[1]]),dim(ugap.core)[1]),file=paste0(prefix,"Core",datestr,".dates.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)


write.table(matrix(as.integer(is.na(ugap.headline)),dim(ugap.headline)[1]),file=paste0(prefix,"Headline",datestr,".yNaN.txt"),sep=" ",row.names=rep("",dim(ugap.headline)[1]),col.names=FALSE,quote=FALSE)
ugap.headline[is.na(ugap.headline)] <- 0
write.table(matrix(sprintf('%30.16E',ugap.headline),dim(ugap.headline)[1]),file=paste0(prefix,"Headline",datestr,".yData.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrix(sprintf('%30.16E',matlab.dates[1:dim(ugap.headline)[1]]),dim(ugap.headline)[1]),file=paste0(prefix,"Headline",datestr,".dates.txt"),sep="",row.names=FALSE,col.names=FALSE,quote=FALSE)
