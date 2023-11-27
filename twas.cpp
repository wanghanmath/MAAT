#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	arma::uvec z_nok_func(arma::uvec z, int k){
		int p=z.size();
		arma::uvec logic=arma::ones<arma::uvec>(p);
		logic(k)=0;
		arma::uvec ids = find(logic == 1);
		arma::uvec z_nok=z.elem(ids);
		return z_nok;
	}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	arma::mat W_nok_func(arma::mat annotation_mat, int k){
		int p=annotation_mat.n_rows;
		arma::uvec logic=arma::ones<arma::uvec>(p);
		logic(k)=0;
		arma::uvec ids = find(logic == 1);
		arma::mat W_nok=annotation_mat.rows(ids);
		return W_nok;
	}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	double g_sim(arma::mat W){
		int m=W.n_cols;
		int n_j=W.n_rows;
		arma::vec g_wl(m);
		double g=0.0, varnew, H;
		double PInew = 3.14159265358979323846;
		for (int l=0;l<m;l++){
			if (n_j==1){
				varnew=1e-2;
                H=0.0;
			}
			if (n_j>1){
				varnew=var(W.col(l));
				if (varnew==0.0){
					varnew=1e-6;
				}
				H=sum(pow((W.col(l)-mean(W.col(l))),2))/(2.0*varnew);
			}
			g_wl(l)=-1.0*n_j*log(2*PInew*varnew)/2.0-H;
			g+=g_wl(l);
		}		
		return g;
	}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	arma::vec update_beta(arma::vec beta, arma::uvec z, arma::mat XHinvX_all, arma::vec sigma2_j, double sigma2_epsilon, arma::mat genotype_center, arma::vec HinvY){
    int p=beta.size();
    double XHinvX,beta_var,xhinvy,beta_mean,xhinvxbeta;
    int z_k;
    arma::vec XHinvY,XHinvXbeta,betanew(p);
    for (int k=0;k<p;k++){
    	XHinvX=XHinvX_all(k,k);
    	z_k=z(k);
        beta_var=sigma2_epsilon/(1.0/sigma2_j(z_k)+XHinvX);
        XHinvY=genotype_center.col(k).t()*HinvY;
        xhinvy=XHinvY(0);
        XHinvXbeta=XHinvX_all.row(k)*beta;
        xhinvxbeta=XHinvXbeta(0)-XHinvX_all(k,k)*beta(k);
        beta_mean=(xhinvy-xhinvxbeta)/(1.0/sigma2_j(z_k)+XHinvX);
        beta(k)=rnorm(1,beta_mean,sqrt(beta_var))(0);
    }
    return beta;    
	}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
    List update_betatest(arma::vec beta, arma::uvec z, arma::mat XHinvX_all, arma::vec sigma2_j, double sigma2_epsilon, arma::mat genotype_center, arma::vec HinvY){
    int p=beta.size();
    double XHinvX,beta_var,xhinvy,beta_mean,xhinvxbeta;
    int z_k;
    arma::vec XHinvY,XHinvXbeta,betanew(p),betamean(p),betavar(p);
    for (int k=0;k<p;k++){
        XHinvX=XHinvX_all(k,k);
        z_k=z(k);
        beta_var=sigma2_epsilon/(1.0/sigma2_j(z_k)+XHinvX);
        XHinvY=genotype_center.col(k).t()*HinvY;
        xhinvy=XHinvY(0);
        XHinvXbeta=XHinvX_all.row(k)*beta;
        xhinvxbeta=XHinvXbeta(0)-XHinvX_all(k,k)*beta(k);
        beta_mean=(xhinvy-xhinvxbeta)/(1.0/sigma2_j(z_k)+XHinvX);
        betanew(k)=R::rnorm(beta_mean,sqrt(beta_var));
        betamean[k]=beta_mean;
        betavar[k]=beta_var;
    }
    return List::create(Named("mean") = betamean, Named("beta_var") = betavar);
    }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	List update_z(arma::uvec z, arma::vec beta, arma::mat annotation_mat, arma::vec sigma2_j, double sigma2_epsilon, double a_k, double b_k){
		int p_size=z.size();
		for (int k=0;k<p_size;k++){
		/*int k=5;*/
		arma::uvec z_nok=z_nok_func(z,k), znok_j_loc, zk_j_loc, zero_loc, nonzero_loc, z_k;
		arma::uvec uniq_z=unique(z_nok);
		arma::mat anno_j, annotation_mat_nok;
		int length_z_nok=uniq_z.size();
        int z_indicator, part2, length_z_k, z_k_length, z_kk_loc_length;
        double part1, part3, sd_j;
        arma::mat g_k_all_matrix(length_z_nok,length_z_nok);
        annotation_mat_nok=W_nok_func(annotation_mat,k);
        arma::vec g_k_all, g_k_tilde_all, g_k_allnew, g_k_tilde_allnew;
        if (length_z_nok>=10){
        	z_indicator=1;
        }
        if (length_z_nok<10){
        	z_indicator=0;
        }
        if (z_indicator==1){
        	arma::vec sigma2_j_new(length_z_nok), g_nok(length_z_nok), k_prob_log(length_z_nok) ;
            arma::vec k_prob_part1(length_z_nok), k_prob(length_z_nok), k_prob_part2(length_z_nok), k_prob_part3(length_z_nok);
            for (int j=0;j<length_z_nok;j++){
            	znok_j_loc=find(z_nok == j);            	
		        anno_j=annotation_mat_nok.rows(znok_j_loc);
		        g_nok(j)=g_sim(anno_j);
            }
            for (int j=0;j<length_z_nok;j++){
            	g_k_all_matrix.col(j)=g_nok;
            }
            for (int j=0;j<length_z_nok;j++){
            	z_k=z;
            	z_k[k]=j;
            	zk_j_loc=find(z_k==j);
            	anno_j=annotation_mat.rows(zk_j_loc);
            	g_k_all_matrix(j,j)=g_sim(anno_j);
            	sd_j=sqrt(sigma2_j(j)*sigma2_epsilon);
            	part1=R::dnorm(beta(k),0,sd_j,FALSE);
                part2=zk_j_loc.size()-1;
                part3=g_k_all_matrix(j,j)-g_nok(j);
                k_prob_part1(j)=part1;
                k_prob_part2(j)=part2;
                k_prob_part3(j)=part3;
                sigma2_j_new(j)=sigma2_j(j);
            }
            zero_loc=find(k_prob_part1==0.0);
            int length_zero_loc=zero_loc.size();
            if (length_zero_loc>0){
            	nonzero_loc=find(k_prob_part1!=0.0);
            	double min_part1=min(k_prob_part1.elem(nonzero_loc));
            	if (min_part1<1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*min_part1;
            	}
            	if (min_part1>=1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*1e-233;
            	}
            }
            double max_part3=max(k_prob_part3);
            if (max_part3<=700.0){
            	for (int j=0;j<length_z_nok;j++){
            		k_prob(j)=k_prob_part1(j)*k_prob_part2(j)*exp(k_prob_part3(j));
            	}
            	k_prob=k_prob/sum(k_prob);
            }
            if (max_part3>700.0){
            	for (int j=0;j<length_z_nok;j++){
            		if (k_prob_part2(j)>0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(k_prob_part2(j))+k_prob_part3(j);
            		}
            		if (k_prob_part2(j)==0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(1e-6)+k_prob_part3(j);
            		}
            	}
            	for (int j=0;j<length_z_nok;j++){
            		if ((max_part3-k_prob_part3(j))>700.0){
            			k_prob(j)=0.0;
            		}
            		if ((max_part3-k_prob_part3(j))<=700.0){
            			arma::vec k_prob_logj=k_prob_log-arma::ones<arma::vec>(length_z_nok)*k_prob_log(j);
            			double k_prob_exp=sum(exp(k_prob_logj));
            			k_prob(j)=1.0/k_prob_exp;
            		}
            	}
            }    
            arma::vec k_prob_sum=arma::zeros<arma::vec>(length_z_nok+1);        
            for (int i=0;i<length_z_nok;i++){
            	k_prob_sum(i+1)=sum(k_prob.subvec(0,i));
            }
            double k_unif=R::runif(0,1);
            arma::uvec k_indicator=arma::zeros<arma::uvec>(length_z_nok);
            for (int i=0;i<length_z_nok;i++){
            	if ((k_unif<=k_prob_sum(i+1))&(k_unif>k_prob_sum(i))){
            		k_indicator(i)=1;
            	}
            }
            arma::uvec z_indicator_loc=find(k_indicator==1);
            z(k)=z_indicator_loc(0);
        }
        if (z_indicator==0){
        	z_k=z;
        	arma::uvec z_kk_loc=find(z_k==z_k(k));
        	z_kk_loc_length=z_kk_loc.size();
        	/*add another indicator*/
        	if (z_kk_loc_length==1){
        	    /*change*/
        		if (z_k(k)!=max(z_k)){
        			arma::vec sigma2_j_replace=sigma2_j;
        			sigma2_j_replace(z_k(k))=sigma2_j(max(z));
        			sigma2_j_replace(max(z))=sigma2_j(z_k(k));
        			sigma2_j=sigma2_j_replace;
        			arma::uvec znew=z;
        			arma::uvec znew_max_loc=find(znew==max(znew));
        			int znew_macloc_length=znew_max_loc.size();
        			znew(znew_max_loc)=arma::ones<arma::uvec>(znew_macloc_length)*z_k(k);
        			z=znew;
        			z_k=z;
        		}
        		arma::uvec uniq_z_k=unique(z_k);
        		length_z_k=uniq_z_k.size();
        		g_k_all=arma::ones<arma::vec>(length_z_k);
        		for (int j=0;j<length_z_k;j++){
        			zk_j_loc=find(z_k==j);
        			anno_j=annotation_mat.rows(zk_j_loc);
            	    g_k_all(j)=g_sim(anno_j);            	    
        		}
        		g_k_tilde_all=g_k_all;
        	}
        	if (z_kk_loc_length>1){
        		arma::uvec uniq_z_k=unique(z_k);
        		length_z_k=uniq_z_k.size();
        		z_k(k)=length_z_k;
        		g_k_all=arma::ones<arma::vec>(length_z_k+1);
        		for (int j=0;j<(length_z_k+1);j++){
        			zk_j_loc=find(z_k==j);
        			anno_j=annotation_mat.rows(zk_j_loc);
            	    g_k_all(j)=g_sim(anno_j);   
        		}
        		g_k_tilde_all=g_k_all;
        	}
        	arma::vec g_k_allnew=g_k_all;
        	g_k_tilde_allnew=g_k_tilde_all;
        	arma::uvec uniq_z_knew=unique(z_k);
        	z_k_length=uniq_z_knew.size();
        	arma::vec k_prob(z_k_length), k_prob_log(z_k_length), k_prob_part1(z_k_length), k_prob_part2(z_k_length), k_prob_part3(z_k_length);
            double IG_sigmaj_one_a=0.5+a_k;
            double IG_sigmaj_one_b=(pow(beta(k),2))/(2.0*sigma2_epsilon)+b_k;
            double gamma_ran=R::rgamma(IG_sigmaj_one_a,IG_sigmaj_one_b);
            double sigma2_new=1.0/gamma_ran;
            arma::vec sigma2_j_new(z_k_length);
            arma::mat g_kj_all_matrix(z_k_length,z_k_length);
            for (int j=0;j<z_k_length;j++){
            	g_kj_all_matrix.col(j)=g_k_allnew;
            }
            for (int j=0;j<z_k_length;j++){
            	unsigned int jnew=(unsigned int)j;
            	if (jnew!=z_k(k)){
            		arma::uvec z_kj=z_k;
            		z_kj(k)=j;
            		arma::uvec z_kj_j_loc=find(z_kj==j);
            		anno_j=annotation_mat.rows(z_kj_j_loc);
            		g_kj_all_matrix(j,j)=g_sim(anno_j);
            		g_kj_all_matrix(z_k(k),j)=0.0;
            	}
            }
            for (int j=0;j<z_k_length;j++){
            	unsigned int jnew=(unsigned int)j;
            	if (jnew!=z_k(k)){
            		arma::uvec z_j_loc=find(z==j);
            	    sd_j=sqrt(sigma2_j(j)*sigma2_epsilon);
            	    part1=R::dnorm(beta(k),0,sd_j,FALSE);
                    part2=z_j_loc.size()-1;
                    part3=g_kj_all_matrix(j,j)-g_k_tilde_allnew(j);
                    k_prob_part1(j)=part1;
                    k_prob_part2(j)=part2;
                    k_prob_part3(j)=part3;
                    sigma2_j_new(j)=sigma2_j(j);
            	}
            	if (jnew==z_k(k)){
            		sd_j=sqrt(sigma2_new*sigma2_epsilon);
            		part1=R::dnorm(beta(k),0,sd_j,FALSE);
            		part2=1;
            		part3=g_k_tilde_allnew(j);
            		k_prob_part1(j)=part1;
            		k_prob_part2(j)=part2;
            		k_prob_part3(j)=part3;
            		sigma2_j_new(j)=sigma2_new;
            	}
            }
            zero_loc=find(k_prob_part1==0.0);
            int length_zero_loc=zero_loc.size();
            if (length_zero_loc>0){
            	nonzero_loc=find(k_prob_part1!=0.0);
            	double min_part1=min(k_prob_part1.elem(nonzero_loc));
            	if (min_part1<1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*min_part1;
            	}
            	if (min_part1>=1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*1e-233;
            	}
            }
            double max_part3=max(k_prob_part3);
            if (max_part3<=700.0){
            	for (int j=0;j<z_k_length;j++){
            		k_prob(j)=k_prob_part1(j)*k_prob_part2(j)*exp(k_prob_part3(j));
            	}
            	k_prob=k_prob/sum(k_prob);
            }
            if (max_part3>700.0){
            	for (int j=0;j<z_k_length;j++){
            		if (k_prob_part2(j)>0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(k_prob_part2(j))+k_prob_part3(j);
            		}
            		if (k_prob_part2(j)==0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(1e-6)+k_prob_part3(j);
            		}
            	}
            	for (int j=0;j<z_k_length;j++){
            		if ((max_part3-k_prob_part3(j))>700.0){
            			k_prob(j)=0.0;
            		}
            		if ((max_part3-k_prob_part3(j))<=700.0){
            			arma::vec k_prob_logj=k_prob_log-arma::ones<arma::vec>(z_k_length+1)*k_prob_log(j);
            			double k_prob_exp=sum(exp(k_prob_logj));
            			k_prob(j)=1.0/k_prob_exp;
            		}
            	}
            }   
            arma::vec k_prob_sum=arma::zeros<arma::vec>(z_k_length+2);
            for (int i=0;i<z_k_length;i++){
            	k_prob_sum(i+1)=sum(k_prob.subvec(0,i));
            }
            double k_unif=R::runif(0,1);
            arma::uvec k_indicator=arma::zeros<arma::uvec>(z_k_length+1);
            for (int i=0;i<z_k_length;i++){
            	if ((k_unif<=k_prob_sum(i+1))&(k_unif>k_prob_sum(i))){
            		k_indicator(i)=1;
            	}
            }
            arma::uvec z_indicator_loc=find(k_indicator==1);
            z(k)=z_indicator_loc(0);
            if (z_indicator_loc(0)==z_k(k)){
            	sigma2_j=sigma2_j_new;
            }
        }
		}
		return List::create(Named("z") = z, Named("sigma2_j") = sigma2_j);
	}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
	List update_part_z(arma::vec loc_loop, arma::uvec z, arma::vec beta, arma::mat annotation_mat, arma::vec sigma2_j, double sigma2_epsilon, double a_k, double b_k){
		int p_size=z.size();
		int loc_size=loc_loop.size();
		for (int loc=0;loc<loc_size;loc++){
        //int loc=0;
		int k=loc_loop(loc);
		arma::uvec z_nok=z_nok_func(z,k), znok_j_loc, zk_j_loc, zero_loc, nonzero_loc, z_k;
		arma::uvec uniq_z=unique(z_nok);
		arma::mat anno_j, annotation_mat_nok;
		int length_z_nok=uniq_z.size();
        int z_indicator, part2, length_z_k, z_k_length, z_kk_loc_length;
        double part1, part3, sd_j;
        arma::mat g_k_all_matrix(length_z_nok,length_z_nok);
        annotation_mat_nok=W_nok_func(annotation_mat,k);
        arma::vec g_k_all, g_k_tilde_all, g_k_allnew, g_k_tilde_allnew;
        if (length_z_nok>=10){
        	z_indicator=1;
        }
        if (length_z_nok<10){
        	z_indicator=0;
        }
        if (z_indicator==1){
        	arma::vec sigma2_j_new(length_z_nok), g_nok(length_z_nok), k_prob_log(length_z_nok) ;
            arma::vec k_prob_part1(length_z_nok), k_prob(length_z_nok), k_prob_part2(length_z_nok), k_prob_part3(length_z_nok);
            for (int j=0;j<length_z_nok;j++){
            	znok_j_loc=find(z_nok == j);            	
		        anno_j=annotation_mat_nok.rows(znok_j_loc);
		        g_nok(j)=g_sim(anno_j);
            }
            for (int j=0;j<length_z_nok;j++){
            	g_k_all_matrix.col(j)=g_nok;
            }
            for (int j=0;j<length_z_nok;j++){
            	z_k=z;
            	z_k[k]=j;
            	zk_j_loc=find(z_k==j);
            	anno_j=annotation_mat.rows(zk_j_loc);
            	g_k_all_matrix(j,j)=g_sim(anno_j);
            	sd_j=sqrt(sigma2_j(j)*sigma2_epsilon);
            	part1=R::dnorm(beta(k),0,sd_j,FALSE);
                part2=zk_j_loc.size()-1;
                part3=g_k_all_matrix(j,j)-g_nok(j);
                k_prob_part1(j)=part1;
                k_prob_part2(j)=part2;
                k_prob_part3(j)=part3;
                sigma2_j_new(j)=sigma2_j(j);
            }
            zero_loc=find(k_prob_part1==0.0);
            int length_zero_loc=zero_loc.size();
            if (length_zero_loc>0){
            	nonzero_loc=find(k_prob_part1!=0.0);
            	double min_part1=min(k_prob_part1.elem(nonzero_loc));
            	if (min_part1<1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*min_part1;
            	}
            	if (min_part1>=1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*1e-233;
            	}
            }
            double max_part3=max(k_prob_part3);
            if (max_part3<=700.0){
            	for (int j=0;j<length_z_nok;j++){
            		k_prob(j)=k_prob_part1(j)*k_prob_part2(j)*exp(k_prob_part3(j));
            	}
            	k_prob=k_prob/sum(k_prob);
            }
            if (max_part3>700.0){
            	for (int j=0;j<length_z_nok;j++){
            		if (k_prob_part2(j)>0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(k_prob_part2(j))+k_prob_part3(j);
            		}
            		if (k_prob_part2(j)==0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(1e-6)+k_prob_part3(j);
            		}
            	}
            	for (int j=0;j<length_z_nok;j++){
            		if ((max_part3-k_prob_part3(j))>700.0){
            			k_prob(j)=0.0;
            		}
            		if ((max_part3-k_prob_part3(j))<=700.0){
            			arma::vec k_prob_logj=k_prob_log-arma::ones<arma::vec>(length_z_nok)*k_prob_log(j);
            			double k_prob_exp=sum(exp(k_prob_logj));
            			k_prob(j)=1.0/k_prob_exp;
            		}
            	}
            }    
            arma::vec k_prob_sum=arma::zeros<arma::vec>(length_z_nok+1);        
            for (int i=0;i<length_z_nok;i++){
            	k_prob_sum(i+1)=sum(k_prob.subvec(0,i));
            }
            double k_unif=R::runif(0,1);
            arma::uvec k_indicator=arma::zeros<arma::uvec>(length_z_nok);
            for (int i=0;i<length_z_nok;i++){
            	if ((k_unif<=k_prob_sum(i+1))&(k_unif>k_prob_sum(i))){
            		k_indicator(i)=1;
            	}
            }
            arma::uvec z_indicator_loc=find(k_indicator==1);
            z(k)=z_indicator_loc(0);
        }
        if (z_indicator==0){
        	z_k=z;
        	arma::uvec z_kk_loc=find(z_k==z_k(k));
        	z_kk_loc_length=z_kk_loc.size();
        	/*add another indicator*/
        	if (z_kk_loc_length==1){
        	    /*change*/
        		if (z_k(k)!=max(z_k)){
        			arma::vec sigma2_j_replace=sigma2_j;
        			sigma2_j_replace(z_k(k))=sigma2_j(max(z));
        			sigma2_j_replace(max(z))=sigma2_j(z_k(k));
        			sigma2_j=sigma2_j_replace;
        			arma::uvec znew=z;
        			arma::uvec znew_max_loc=find(znew==max(znew));
        			int znew_macloc_length=znew_max_loc.size();
        			znew(znew_max_loc)=arma::ones<arma::uvec>(znew_macloc_length)*z_k(k);
        			z=znew;
        			z_k=z;
        		}
        		arma::uvec uniq_z_k=unique(z_k);
        		length_z_k=uniq_z_k.size();
        		g_k_all=arma::ones<arma::vec>(length_z_k);
        		for (int j=0;j<length_z_k;j++){
        			zk_j_loc=find(z_k==j);
        			anno_j=annotation_mat.rows(zk_j_loc);
            	    g_k_all(j)=g_sim(anno_j);            	    
        		}
        		g_k_tilde_all=g_k_all;
        	}
        	if (z_kk_loc_length>1){
        		arma::uvec uniq_z_k=unique(z_k);
        		length_z_k=uniq_z_k.size();
        		z_k(k)=length_z_k;
        		g_k_all=arma::ones<arma::vec>(length_z_k+1);
        		for (int j=0;j<(length_z_k+1);j++){
        			zk_j_loc=find(z_k==j);
        			anno_j=annotation_mat.rows(zk_j_loc);
            	    g_k_all(j)=g_sim(anno_j);   
        		}
        		g_k_tilde_all=g_k_all;
        	}
        	arma::vec g_k_allnew=g_k_all;
        	g_k_tilde_allnew=g_k_tilde_all;
        	arma::uvec uniq_z_knew=unique(z_k);
        	z_k_length=uniq_z_knew.size();
        	arma::vec k_prob(z_k_length), k_prob_log(z_k_length), k_prob_part1(z_k_length), k_prob_part2(z_k_length), k_prob_part3(z_k_length);
            double IG_sigmaj_one_a=0.5+a_k;
            double IG_sigmaj_one_b=(pow(beta(k),2))/(2.0*sigma2_epsilon)+b_k;
            double gamma_ran=R::rgamma(IG_sigmaj_one_a,IG_sigmaj_one_b);
            double sigma2_new=1.0/gamma_ran;
            arma::vec sigma2_j_new(z_k_length);
            arma::mat g_kj_all_matrix(z_k_length,z_k_length);
            for (int j=0;j<z_k_length;j++){
            	g_kj_all_matrix.col(j)=g_k_allnew;
            }
            for (int j=0;j<z_k_length;j++){
            	unsigned int jnew=(unsigned int)j;
            	if (jnew!=z_k(k)){
            		arma::uvec z_kj=z_k;
            		z_kj(k)=j;
            		arma::uvec z_kj_j_loc=find(z_kj==j);
            		anno_j=annotation_mat.rows(z_kj_j_loc);
            		g_kj_all_matrix(j,j)=g_sim(anno_j);
            		g_kj_all_matrix(z_k(k),j)=0.0;
            	}
            }
            for (int j=0;j<z_k_length;j++){
            	unsigned int jnew=(unsigned int)j;
            	if (jnew!=z_k(k)){
            		arma::uvec z_j_loc=find(z==j);
            	    sd_j=sqrt(sigma2_j(j)*sigma2_epsilon);
            	    part1=R::dnorm(beta(k),0,sd_j,FALSE);
                    part2=z_j_loc.size()-1;
                    part3=g_kj_all_matrix(j,j)-g_k_tilde_allnew(j);
                    k_prob_part1(j)=part1;
                    k_prob_part2(j)=part2;
                    k_prob_part3(j)=part3;
                    sigma2_j_new(j)=sigma2_j(j);
            	}
            	if (jnew==z_k(k)){
            		sd_j=sqrt(sigma2_new*sigma2_epsilon);
            		part1=R::dnorm(beta(k),0,sd_j,FALSE);
            		part2=1;
            		part3=g_k_tilde_allnew(j);
            		k_prob_part1(j)=part1;
            		k_prob_part2(j)=part2;
            		k_prob_part3(j)=part3;
            		sigma2_j_new(j)=sigma2_new;
            	}
            }
            zero_loc=find(k_prob_part1==0.0);
            int length_zero_loc=zero_loc.size();
            if (length_zero_loc>0){
            	nonzero_loc=find(k_prob_part1!=0.0);
            	double min_part1=min(k_prob_part1.elem(nonzero_loc));
            	if (min_part1<1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*min_part1;
            	}
            	if (min_part1>=1e-233){
            		k_prob_part1(zero_loc)=arma::ones<arma::vec>(length_zero_loc)*1e-233;
            	}
            }
            double max_part3=max(k_prob_part3);
            if (max_part3<=700.0){
            	for (int j=0;j<z_k_length;j++){
            		k_prob(j)=k_prob_part1(j)*k_prob_part2(j)*exp(k_prob_part3(j));
            	}
            	k_prob=k_prob/sum(k_prob);
            }
            if (max_part3>700.0){
            	for (int j=0;j<z_k_length;j++){
            		if (k_prob_part2(j)>0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(k_prob_part2(j))+k_prob_part3(j);
            		}
            		if (k_prob_part2(j)==0){
            			k_prob_log(j)=log(k_prob_part1(j))+log(1e-6)+k_prob_part3(j);
            		}
            	}
            	for (int j=0;j<z_k_length;j++){
            		if ((max_part3-k_prob_part3(j))>700.0){
            			k_prob(j)=0.0;
            		}
            		if ((max_part3-k_prob_part3(j))<=700.0){
            			arma::vec k_prob_logj=k_prob_log-arma::ones<arma::vec>(z_k_length+1)*k_prob_log(j);
            			double k_prob_exp=sum(exp(k_prob_logj));
            			k_prob(j)=1.0/k_prob_exp;
            		}
            	}
            }   
            arma::vec k_prob_sum=arma::zeros<arma::vec>(z_k_length+2);
            for (int i=0;i<z_k_length;i++){
            	k_prob_sum(i+1)=sum(k_prob.subvec(0,i));
            }
            double k_unif=R::runif(0,1);
            arma::uvec k_indicator=arma::zeros<arma::uvec>(z_k_length+1);
            for (int i=0;i<z_k_length;i++){
            	if ((k_unif<=k_prob_sum(i+1))&(k_unif>k_prob_sum(i))){
            		k_indicator(i)=1;
            	}
            }
            arma::uvec z_indicator_loc=find(k_indicator==1);
            z(k)=z_indicator_loc(0);
            if (z_indicator_loc(0)==z_k(k)){
            	sigma2_j=sigma2_j_new;
            }
        }
		}
		return List::create(Named("z") = z, Named("sigma2_j") = sigma2_j);
	}













