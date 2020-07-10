# Data Preprocessing ------------------------------------------------------

load('BRCA/BRCA_Data_COBS.RData')

dataY <- cbind(scale(Expr_large,center=T,scale=T),
            scale(Meth_large,center=T,scale=T),
            scale(CNV_large,center=T,scale=T))
dataX <- BRCA_Tumor  

# Retrieve the dimensions of the data matrices
n=dim(dataY)[1]; p=dim(dataY)[2]; q=dim(dataX)[2]
gp=ceiling(1:p/200) #group index of Y, multi-block 
pvec <- vector()  
for (i in 1:length(unique(gp))) {
  pvec[i] <- sum( gp == unique(gp)[i] ) 
}
gpb = rep(1,q)
qvec <- vector()  
for (i in 1:length(unique(gpb))) {
  qvec[i] <- sum( gpb == unique(gpb)[i] ) 
}

cat('BRCA data \n')
cat("Sample size: ", n , "\n")
cat("Number of variables", p, ", with", length(unique(gp)), "blocks of sizes", pvec,"\n")
cat("Number of covariates", q, ", with", length(unique(gpb)), "blocks of sizes", qvec,"\n")

# Check that the data is already centered; 
cat("Primary data (Y) is centered if", max(abs(apply(dataY,2,mean))), "is small","\n")

 
 
# rownames(dataY) <- NULL; rownames(dataX)=NULL; colnames(dataY)=NULL;
 
# X = dataX (dim n x q)
# Y = dataY (dim n x p)
# Y_index = gp (group indice given by numeric values)
# X_index = null (no group exists)

rm("CNV_large","Expr_large","Meth_large","BRCA_Tumor","i")
