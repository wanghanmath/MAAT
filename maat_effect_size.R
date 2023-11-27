library("Rcpp")
library("optparse")
library("mvtnorm")
library("mclust")
library("glmnet")
library("MASS")
sourceCpp("twas.cpp")
args <- commandArgs(TRUE)
genotype_path<-paste0(args[1],args[2])
expression_path<-paste0(args[1],args[3])
annotation_path<-paste0(args[1],args[4])
genotype<-read.csv(genotype_path,header=T,stringsAsFactors=F)
expression<-read.csv(expression_path,header=T,stringsAsFactors=F)
expression<-c(expression[,1])
annotation_matrix<-read.csv(annotation_path,header=T,stringsAsFactors=F)
sparsity_quant<-args[5]
setwd(args[1])
genotype_center<-scale(genotype,scale=T)
expression<-scale(expression,scale=T)
#set initial values for z based on gaussian mixture model in mclust#
anno_mclust<-Mclust(annotation_matrix[,1:7])
reorder<-rank(anno_mclust$parameters$mean[1,])
class<-rep(0,dim(genotype)[2])
for (i in 1:length(reorder)){
    classification_uniq<-sort(unique(anno_mclust$classification))
    uniq_i<-classification_uniq[i]
	class[which(anno_mclust$classification==uniq_i)]<-reorder[i]
}
max_class<-max(class)
class[which(class==0)]<-max_class+1
h2<-0.2
sigma02<-0.3
n<-dim(genotype)[1]
p<-dim(genotype)[2]
#set prior parameters#
a_k<-b_k<-0.1
a_epsilon<-b_epsilon<-0.1
a_beta<-2
b_beta<-8
class_num<-length(unique(class))
f_h<-function(h){
	residual<-as.matrix(expression-genotype_center%*%as.matrix(beta))
	part1<-log(1/sqrt(det(H)))
	part2<-(-1/(2*sigma2_epsilon))*t(residual)%*%ginv(H)%*%residual
	part3<-log(((1-h)/h)^(a_k+1))
	part4=-b_k*(1-h)/h
	part5<-log(1/(1-h)^2)
	f_h<-part1+part2+part3+part4+part5
	return(f_h)
}	
#set initial values#
beta_loop<-matrix(ncol=200,nrow=dim(genotype)[2])
z_loop<-matrix(ncol=200,nrow=dim(genotype)[2])
sigma2_zero_loop<-sigma2_epsilon_loop<-rep(0,200)
sigma2_j_loop<-matrix(nrow=dim(genotype)[2],ncol=200)
lasso<-glmnet(genotype_center,expression)
beta<-c(lasso$beta[,dim(lasso$beta)[2]])
zscore_snp<-rep(0,dim(genotype)[2])
for (i in 1:length(zscore_snp)){
    aa<-lm(expression~genotype_center[,i])
    aa_sum<-summary(aa)
    coef<-aa_sum$coefficients[2,1]
    sd<-aa_sum$coefficients[2,2]
    zscore_snp[i]<-coef/sd
}
update_loc<-c(which(beta!=0))
if (length(update_loc)<p){
    order_zscore<-order(abs(zscore_snp),decreasing=T)
    order_reserved<-setdiff(order_zscore,update_loc)
    reserved_length<-p-length(update_loc)
    update_locnew<-order_reserved[1:reserved_length]
    update_loc<-c(update_loc,update_locnew)
}
loc_loop<-as.numeric(sort(update_loc))
loc_loop<-loc_loop-1
z<-class
z<-z-1
sigma2_zero<-0.1
h<-sigma2_zero/(1+sigma2_zero)
alpha<-rep(0,dim(annotation_matrix)[2])
for (i in 1:length(alpha)){
	alpha[i]<-var(annotation_matrix[,i])
}
residual<-as.matrix(expression-genotype_center%*%as.matrix(beta))
sigma2_epsilon<-t(residual)%*%residual/(n-1)  #similar to DPR#
sigma2_j<-rep(0,length(unique(class)))
for (j in 1:length(unique(class))){
	IG_sigmaj_a<-length(which(class==j))/2+a_k
	beta_square<-beta*beta
	IG_sigmaj_b<-sum(beta_square[which(class==j)])/(2*sigma2_epsilon)+b_k
    sigma2_j[j]<-1/rgamma(1,IG_sigmaj_a,IG_sigmaj_b)
}
beta_loop[,1]<-beta
z_loop[,1]<-z
sigma2_zero_loop[1]<-sigma2_zero
sigma2_epsilon_loop[1]<-sigma2_epsilon
sigma2_j_loop[1:length(sigma2_j),1]<-c(sigma2_j)
annotation_matrix<-as.matrix(annotation_matrix)
for (iter in 2:200){
    #update beta#
    H<-diag(n)+sigma2_zero*genotype_center%*%t(genotype_center)/p
    Hinv<-ginv(H)
    HinvY<-Hinv%*%as.matrix(expression)
    XHinvX_all<-t(genotype_center)%*%Hinv%*%genotype_center
    beta<-update_beta(beta, z, XHinvX_all, sigma2_j, sigma2_epsilon, genotype_center, HinvY)
    beta_loop[,iter]<-c(beta)
    #update z#
    if (iter%%5==0){
        z_list<-update_z(z,beta,annotation_matrix,sigma2_j,sigma2_epsilon,a_k,b_k)
        z<-z_list$z
        sigma2_j<-z_list$sigma2_j
    }
    if (iter%%5!=0){
        z_list<-update_part_z(loc_loop,z,beta,annotation_matrix,sigma2_j,sigma2_epsilon,a_k,b_k)
        z<-z_list$z
        sigma2_j<-z_list$sigma2_j
    }   
    z_loop[,iter]<-c(z)
    #update sigma2_j#
    sigma2_j<-rep(0,length(unique(z)))
    for (j in 1:length(unique(z))){
        IG_sigmaj_a<-length(which(z==(j-1)))/2+a_k
        beta_square<-beta*beta
        IG_sigmaj_b<-sum(beta_square[which(z==(j-1))])/(2*sigma2_epsilon)+b_k
        sigma2_j[j]<-1/rgamma(1,IG_sigmaj_a,IG_sigmaj_b)
    }
    sigma2_j_loop[c(1:length(sigma2_j)),iter]<-c(sigma2_j)
    #update sigma2_epsilon#
    IG_sigma_epsilon_a<-n/2+a_epsilon+p/2
    Y_Xbeta<-as.matrix(expression-genotype_center%*%as.matrix(beta))
    IG_sigma_epsilon_b1<-0.5*t(Y_Xbeta)%*%Hinv%*%Y_Xbeta
    IG_sigma_epsilon_b2<-b_epsilon
    IG_sigma_epsilon_b3<-0
    for (j in 1:length(unique(z))){
        IG_sigma_epsilon_b3<-IG_sigma_epsilon_b3+sum(beta_square[which(z==(j-1))])/(2*sigma2_j[j])
    }
    IG_sigma_epsilon_b<-IG_sigma_epsilon_b1+IG_sigma_epsilon_b2+IG_sigma_epsilon_b3
    sigma2_epsilon<-1/rgamma(1,IG_sigma_epsilon_a,IG_sigma_epsilon_b)
    sigma2_epsilon_loop[iter]<-sigma2_epsilon
    #update sigma2_zero#
    h2_new<-rbeta(1,a_beta,b_beta)
    mh_ratio1<-exp(f_h(h2_new)-f_h(h))
    mh_ratio2<-dbeta(h,a_beta,b_beta)/dbeta(h2_new,a_beta,b_beta)
    mh_ratio<-mh_ratio1*mh_ratio2
    mh_ratio_min<-min(mh_ratio,1)
    h_unif<-runif(1,0,1)
    if (h_unif<mh_ratio_min){
        h<-h2_new
    }
    sigma2_zero<-h/(1-h)
    sigma2_zero_loop[iter]<-sigma2_zero
}
#post-processing beta to a given sparsity level#
beta_maat<-apply(beta_loop[,101:200],1,mean)
indicator_matrix<-matrix(0,ncol=100,nrow=p)
for (m in 101:200){
    sigma2_j_order<-order(sigma2_j_loop[1:10,m],decreasing=F)-1
    indicator<-c()
    for (kk in 1:2){
        indicator_set<-which(z_loop[,m]==sigma2_j_order[kk])
        indicator<-c(indicator,indicator_set)
    }
    indicator_matrix[c(indicator),(m-100)]<-1
}
beta_indicator<-c(apply(indicator_matrix,1,sum))
write.csv(beta_loop,file=paste0(args[1],"beta_loop.csv"),row.names=F,quote=F)
write.csv(z_loop,file=paste0(args[1],"z_loop.csv"),row.names=F,quote=F)
write.csv(sigma2_j_loop,file=paste0(args[1],"sigma2_j_loop.csv"),row.names=F,quote=F)
write.csv(sigma2_epsilon_loop,file=paste0(args[1],"sigma2_epsilon_loop.csv"),row.names=F,quote=F)
write.csv(sigma2_zero_loop,file=paste0(args[1],"sigma2_zero_loop.csv"),row.names=F,quote=F)
drop_indicator<-which(beta_indicator>=quantile(beta_indicator,sparsity_quant))
beta_maat[drop_indicator]<-0
write.csv(beta_maat,file=paste0(args[1],"beta_maat.csv"))

