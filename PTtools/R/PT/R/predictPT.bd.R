.packageName <- 'PT'

# finiteN prediction of anisotropic undersampling for
# "Sparsity/Undersampling Tradeoffs in Anisotropic Undersampling..."
# by Monajemi and Donoho(I&I 2018)


fprintf <- function(...){
cat(sprintf(...))
}

predictPT.bd <- function(delta, field, N, B=NULL, qRegular=TRUE, so=FALSE){
    
    if(is.null(B)){
        B = N;
    }
    

Eps0  = predictPT(delta, field);

gamma0 = sqrt( 2 * log(B) / N );

res = eta.bd(delta,field);
eta.delta    = res$eta;
alpha   	 = res$alpha;
beta    	 = res$beta;


   
   # compute fractional capacity deficit  R.w = (eps* - epsN)/eps*
   if(!so){
       R.W    =  gamma0  * (alpha * eta.delta);                            	 # reg, first order
   }else if(so & field=="Bnd"){
       R.W    = gamma0 *sqrt( (alpha* eta.delta)^2 + (beta * gamma0/Eps0)^2 );  #excat for Bnd
   }else if(so){
       R.W    =   eta.delta *  (alpha * gamma0 + beta * gamma0^2);       # reg, second order
   }
   
   # compute finiteN epsilon
   
   
    EpsN  = Eps0 * (1-R.W);
    
    EpsN[EpsN<0] = 0;

    # for irregular case, add extra offset due to irregularity
   if(!qRegular)
      EpsN = EpsN - gamma0 * sqrt(EpsN);                          # lower Bound
     #EpsN = EpsN - 0.3678794/sqrt(alpha) * gamma0 * sqrt(EpsN);  # prediction .
   end



return(EpsN)
}




