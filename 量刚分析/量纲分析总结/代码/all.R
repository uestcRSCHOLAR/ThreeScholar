#   system.time(source("all.R",encoding="utf-8"))
#   database <- read.csv('DimensionalDatabase.csv', header=T, as.is=T)  #量纲数据库
#   data <- read.csv('DataTest_all.csv', header=T, as.is=T)   #数据from java

  setwd("C:/Users/CISDI/Desktop/bug/量纲分析/测试")
  data <- read.csv('量纲分析.csv', header=T, as.is=T)   #数据from java
  labelsRAW <- c("",  "mm",  "t", "mm", "t",  "mm", "mm", "C", "C", "C", "m/min", "m/s", "m/s", "kpa", "kpa", "kpa", "C", "C", "C") #java传过来的数据
  transformv <- c(1, 0.001, 1000, 0.001, 1000, 0.001, 0.001, 1, 1, 1, 0.016667,1,1,1000,1000,1000,1,1,1) #java传过来的数据
  DatabaseUnitLabel <-c("","L_DB","M_DB","L_DB","M_DB","L_DB","L_DB","Seta_DB","Seta_DB","Seta_DB","Vs_DB","Vs_DB","Vs_DB","Pr_DB","Pr_DB","Pr_DB","Seta_DB","Seta_DB","Seta_DB" )
#--------------------------------------------------------------------------------量纲特征库
 databaseEle <- colnames(read.csv('electricitydatabase.csv', header=T, as.is=T))      #电学量纲特征库列名称
 databaseTher <- colnames(read.csv('thermodynamicsdatabase.csv', header=T, as.is=T))  #热学量纲特征库列名称
 databasepower <- colnames(read.csv('powerdatabase.csv', header=T, as.is=T))          #力学量纲特征库列名称
#--------------------------------------------------------------------------------量纲特征库
#   labelsRAW <- c("", "cm", "km", "t", "kg", "m/min", "m/s", "kpa", "kpa", "pa", "J", "J", "kJ") #java传过来的数据
#   transformv <- c(1, 0.01, 1000, 1000, 1, 0.016667, 1, 1000, 1000, 0.001, 1, 1, 1000) #java传过来的数据
#   DatabaseUnitLabel <-c("", "L_DB", "L_DB","M_DB", "M_DB", "Vs_DB", "Vs_DB",  "Pr_DB", "Pr_DB", "Pr_DB", "En_DB", "En_DB", "En_DB")
 increasedata <- main(data,  labelsRAW, transformv ,DatabaseUnitLabel)

#--------------------------------------------------------------------数据净化器
 naresult <- apply(increasedata, 2, function(x){ sum(is.na(x)) })
 if ( length(which(naresult !=0)) > 0 ){
   increasedata <- increasedata[,-which(naresult !=0)]
}

 nanresult <- apply(increasedata, 2, function(x){ sum(is.nan(x)) })
 if ( length(which(nanresult !=0)) > 0 ){
    increasedata <- increasedata[,-which(nanresult !=0)]
}

 Infresult <- apply(increasedata, 1, function(x){ sum(is.infinite(x))})
 if ( length(which(Infresult !=0)) > 0 ){
    increasedata <- increasedata[,-which(Infresult !=0)]
}
#--------------------------------------------------------------------数据净化器



 write.table(increasedata ,file='C:/Users/CISDI/Desktop/bug/量纲分析/测试/resultall.csv',row.names=F,sep=",",quote =TRUE )

main <- function( data,  labelsRAW, transformv ,DatabaseUnitLabel) {
#--------------------------------------------------------------------------------根据物理量自动匹配库
  uniquelabel <- unique(DatabaseUnitLabel)
  if (length(intersect(DatabaseUnitLabel,databaseEle))>=2 | length(intersect(DatabaseUnitLabel,databaseTher))>=2) {

    if (length(intersect(DatabaseUnitLabel,databaseEle))>=2 & length(intersect(DatabaseUnitLabel,databaseTher))>=2) {
       database<- read.csv('comlexdtabase.csv', header=T, as.is=T)  #读取综合库
       labels <-intersect(DatabaseUnitLabel, colnames(database))    #用于求解量刚方程组
       if(length(uniquelabel)-length(which(uniquelabel=="")) < 6){
        stop("error:physical symbolic is less than five")
      }
    } else {
      database <- read.csv('comlexdtabase1.csv', header=T, as.is=T) #读取综合库1
      labels <-intersect(DatabaseUnitLabel, colnames(database))      #用于求解量刚方程组
      if(length(uniquelabel)-length(which(uniquelabel=="")) < 5){
        stop("error:physical symbolic is less than five")
      }
    }
  } else  {
      database <-read.csv('Forcedatabase.csv', header=T, as.is=T)    #读取力学库
      labels <-intersect(DatabaseUnitLabel, colnames(database))      #用于求解量刚方程组

      if (length(uniquelabel)-length(which(uniquelabel=="")) < 4){
          stop("error:physical symbolic is less than four")
      }

  }
  #--------------------------------------------------------------------------------根据物理量自动匹配库

  #save.image("test.RData")
  #load("test.RData")
  #提取有单位的数据并进行进制的转换得到处理后的data----------------------------------------------------
  rawdataNames <-colnames(data)          #需要增维的数据列名称
  sondatau <- which(labelsRAW != "")     #有单位的位置
  DatabaseUnitLabel <- DatabaseUnitLabel[sondatau] #得到新数据得到物理量符号
  sonData <- as.matrix(data[, sondatau]) #有单位的数据
  sondataT <- transformv[sondatau]       #转化需要的进度
  diagsonData <- diag(sondataT)  #制造对角矩阵
  data <-  sonData %*% diagsonData
  colnames(data) <- rawdataNames[sondatau]
  #提取有单位的数据并进行进制的转换得到处理后的data----------------------------------------------------
  #write.table(data ,file='C:/Users/CISDI/Desktop/bug/量纲分析/测试/clt.csv',row.names=F,sep=",",quote =TRUE )

  dataA <- as.matrix(subset(database, select=labels)) #提取解方程组需要的量纲数据
  allover1 <- main2(dataA=dataA, data=data,labels=labels, DatabaseUnitLabel=DatabaseUnitLabel,nega=1)     #dataA 用于解方程组，data用于增维
  allover2 <- main2(dataA=dataA, data=data,labels=labels, DatabaseUnitLabel=DatabaseUnitLabel,nega=-1)     #dataA 用于解方程组，data用于增维
  allover3 <- cbind( allover1, allover2)
  allover4 <- Increasedata2( data=data, uniquelabel=uniquelabel, DatabaseUnitLabel=DatabaseUnitLabel) #同一个单位内部的循环
  allover <- cbind( allover3, allover4)

  return(allover)

   #write.table(allover ,file='C:/Users/CISDI/Desktop/bug/量纲分析/测试/Increase.csv',row.names=F,sep=",",quote =TRUE )
   #本函数主要功能：用于适合哪个量纲库的条件判断
   #主函数
  }

cltsolve <- function(dataA, tol) {
  m <- dim(dataA)[1] #dataA矩阵对应的行数
  n <- dim(dataA)[2] #dataA矩阵对应的列数
  i <- 1
  j <- 1
  jb <- c()

  #这里while是做初等变换为下三角
  while (i <= m & j <= n) {


    Media <- abs(as.matrix(dataA[i:m,j]))
    p <- max( Media)       #取出列j对应的最大值
    k <- which(Media == p) #列j最大值对应的位置
    k <- k+i-1  #为后续更新取最大值位置做准备

    if (p <= tol) {
      dataA[i:m,j]  <- matrix(0,m-i+1,1)  #如果小于给予的精度就设置为0
      j <- j + 1
    }
    else {
      jb <- c(jb, j) #储存基向量的位置
      dataA[c(i, k), j:n] <- dataA[c(k,i), j:n] #换行最大值所在第k行和第i行交换
      dataA[i,j:n] <- dataA[i,j:n]/dataA[i,j]
      #------------------制造循环变量
      c1 <- 1:(i-1)
      if (0 %in% c1){
        c1<- NULL
      }
      c2 <- (i+1):m
      if (i+1 > m) {
        c2 <- NULL
      }
      CC <- c(c1,c2)
      #------------------制造循环变量

      for (k in  CC) {

        dataA[k,j:n] <- dataA[k,j:n] - dataA[k,j]*dataA[i,j:n]

      }

      i <- i + 1
      j <- j + 1
    }

  }

  pivcol <- jb       #基向量的标签
  R <- dataA
  r <- length(pivcol) #秩的个数
  nopiv <- 1:n        #列数
  nopiv <-nopiv[-pivcol]  #选取非基向量对应的位置

  Z <- matrix(0,n,n-r)

  if (n > r) {

    ZZ <- diag(1,n-r)   #生成对角线全为1的对角矩阵
    Z[nopiv, ] <- ZZ
  }

  if (r > 0) {
    Z[pivcol, ] <- -R[1:r, nopiv]
  }

  return(Z)

  #本程序用于线性方程组的求解（线性变化求下三角和求解方程）
}

loop <- function(datagroup) {
  gr <- list()
  i=1
  j=1
  if (length(datagroup)%%2 == 0) {
    while (i <= (length(datagroup)-1)) {
      gr[[j]] <- c(i,i+1)
      i=i+2
      j=j+1

    }} else {
      while (i <= (length(datagroup)-2)) {
        gr[[j]] <- c(i,i+1)
        i=i+2
        j=j+1
      }
      gr[(length(datagroup)+1) / 2] <- length(datagroup)
    }

  result <-list()
  if (length(gr[[length(gr)]])==2){
    for( k in 1:length(gr)) {
      blocka <- gr[[k]][1]
      blockb <- gr[[k]][2]
      block_a <- as.data.frame(datagroup[[blocka]])
      block_a_names <- colnames(block_a)
      block_b <- datagroup[[blockb]]
      block_b_names <- colnames(block_b)
      datalist <- lapply(block_a,function(x){x*block_b})
      datapart <- do.call("cbind", datalist)  #数据
      datalist_names<-c()
      for (i in 1:length(block_a_names)){
        medianames <- lapply(block_b_names,function(x){paste(block_a_names[i], '*', x)}) #先将数据的列名称*到一起
        medianames <- do.call("cbind",medianames)
        datalist_names[[i]] <- medianames   #数据列名称
      }
      datapart_names <- do.call("cbind", datalist_names)
      colnames(datapart) <- datapart_names[1,]
      result[[k]]<- datapart
    }} else {
      for (k in 1:(length(gr)-1)) {
        blocka <- gr[[k]][1]
        blockb <- gr[[k]][2]
        block_a <- as.data.frame(datagroup[[blocka]])
        block_a_names <- colnames(block_a)
        block_b <- datagroup[[blockb]]
        block_b_names <- colnames(block_b)
        datalist <- lapply(block_a,function(x){x*block_b})
        datapart <- do.call("cbind", datalist)  #数据
        datalist_names<-c()
        for (i in 1:length(block_a_names)){
          medianames <- lapply(block_b_names,function(x){paste(block_a_names[i], '*', x)})
          medianames <-do.call("cbind",medianames)
          datalist_names[[i]] <- medianames   #数据列名称
        }
        datapart_names <- do.call("cbind", datalist_names)
        colnames(datapart) <- datapart_names[1,]
        result[[k]]<- datapart
      }
      result[[length(gr)]] <- datagroup[[length(datagroup)]]
    }

  return(result)

#  本函数用于基础解析一个解，来对数据进行二分法组合。
#  列名字的组合
}

Increasedata <- function(number, Z, labels, data, DatabaseUnitLabel) {
  Z<-as.matrix(Z)
  math <- Z[which(Z[,number]!=0), number]    #基础解析矩阵中的第number列
  dataLabels <- labels[which(Z[,number]!=0)] #第number列里面非0数值对应的物理量名称
  dataindex <- lapply(dataLabels,grep,DatabaseUnitLabel)           #数据分组依据的物理量名称
  datagroup <- lapply(dataindex,function(x){datap=data[,x]})       #数据按照物理量名称分组
  options(digits = 5) #数据的小数点
  #数据按照基础解析的值进行处理----------------------------
  mathdata <-list()
  for (j in 1:length(math)) {
    media <- datagroup[[j]]
    medialabel <- math[j]
    mathdata[[j]]<- media^medialabel
  }
  #数据按照基础解析的值进行处理----------------------------
  result_all <- loop(mathdata)
  while (length(result_all)!=1) {
    datagroup <-  result_all
    result_all <- loop(datagroup=datagroup)
  }
  result_all <- result_all[[1]]
  return(result_all)
  #本函数是针对不同单位组合
  #调用loop函数
}

Increasedata2 <- function( data, uniquelabel, DatabaseUnitLabel) {

  if  (length(which(uniquelabel==""))==0){
    dataLabels<- uniquelabel
  }else {
    dataLabels <- uniquelabel[-which(uniquelabel=="")]
  }
  dataindex <- lapply(dataLabels, grep, DatabaseUnitLabel)
  datagroup <- lapply(dataindex, function(x){datap=data[,x]})       #数据按照物理量名称分组
  singal <- list()
  datalist_names <- list()
  for (ii in 1:length(datagroup))  {
    #--------------------------------------------------------------获取数据
    block_a <- as.data.frame(datagroup[[ii]])
    block_a_names <- colnames(block_a)
    n <- ncol( block_a )#列数
    selectCol <- seq(from=1, to=(1+(n-1)*(n+1)), by=(n+1))
    medianames <- lapply(block_a_names,function(x){paste(x, '/',block_a_names)}) #先将数据的列名称/到一起
    matrixMedia <- lapply( medianames,function(x){as.matrix(x)})
    blockMedia <- t(do.call("rbind", matrixMedia))
    cbindmedia  <- gsub(pattern = "[[:blank:]]",replacement = "",blockMedia) #去掉空格
    buttonDatanames <-  cbindmedia[, -selectCol]
    datalist <- lapply(block_a,function(x){x/block_a})
    buttonData <- do.call("cbind", datalist)  #数据
    buttonData <- buttonData[, -selectCol]
    colnames(buttonData) <-  buttonDatanames
    singal[[ii]] <- do.call("cbind", buttonData)  #数据
    datalist <-list()
    medianames <- list()
    #---------------------------------------------------------------获取数据与名字
  }

  interData <- do.call("cbind", singal)

  return( interData )
 #本函数是针对所有单位相同的做组合

 }

main2 <- function (dataA, data, labels, DatabaseUnitLabel, nega) {

  Z <- cltsolve(dataA=dataA, tol=0.00001)  #解方程得到基础解析
  Z <- apply(Z,2,function(x){x*nega})

  count <- ncol(Z)       #基础解析中解的个数
  dataincrease <- list() #储存增维结果
  for (number in 1:count) {
    mm <- Increasedata(number=number, Z=Z, labels=labels, data=data, DatabaseUnitLabel=DatabaseUnitLabel)  #利用函数Increasedatadui基础解析得到的单个解，对应到数据去增维
    dataincrease[[number]] <- mm

    #提取相乘过后的列名称，并对其进行最后的处理得到增维数据的列名称--------------------------
    dataincrease_names <- colnames(dataincrease[[number]]) #提取数据列名称
    math <- Z[which(Z[,number]!=0), number]                #基础解析第number列的值
    printnames <- strsplit(dataincrease_names, "[*]")      #拆分字符串中间有*号的都分开
    pasteStep1=list()
    for (i in 1:length(printnames)){
      pasteStep2=c()
      for  (j in 1:length(math)){
        pasteStep2[j] <- paste(printnames[[i]][j],"^","(",math[j],")")
      }
      pasteStep1[[i]]<-pasteStep2
    }
    pasteStep3 <- list()
    for (h in 1:length(pasteStep1)) {
      A <- pasteStep1[[h]][1]
      for (k in 1:(length(pasteStep1[[h]])-1)){
        A <- paste(A,"*",pasteStep1[[h]][k+1])
      }
      A <- gsub(pattern = "[[:blank:]]",replacement = "", A)
      pasteStep3[[h]] <- A
    }
    pasteStep3_name <- do.call("cbind",pasteStep3)
    colnames(dataincrease[[number]]) <- pasteStep3_name[1,]
  }
  #提取相乘过后的列名称，并对其进行最后的处理得到增维数据的列名称--------------------------

  allover <- do.call("cbind",dataincrease) #增维得到的列表数据转化为矩阵
  return(allover)
  #本函数对全部数据进行和并
  #处理列名
  #调用Increasedata
}























