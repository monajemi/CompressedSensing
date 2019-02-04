.packageName <- 'PT'
# eta function for block diagonal 
# phase transition predcition 
# 
# (eps* - eps.N) / eps^* = alpha * eta(fld,delta) * gamma(N,B) + O(gamma^2) 
# slope = eta * eps^* ;
# eps.N = eps* - slope * (alpha * gamma + beta^gamma^2)
# "Sparisty/Undersamping tradeoffs in Anisotropic Undersampling"
# by Monajemi and Donoho (I&I 2018)


eta.bd <- function(deltas,field){
	
	e.deltas = predictPT(deltas, field);
	r.deltas = e.deltas/deltas;
    
	if(field=="Pos"){
	    f.deltas    = (1/deltas - r.deltas)^.5;
	    alpha = 1.0;
	    beta  = -1/3;
	}else if (field=="R"){
	    f.deltas = sqrt(1/deltas);
	    alpha = 1.0;
	    beta  = -1/2
	}else if (field=="C"){
	    f.deltas = sqrt(1/deltas);
	    alpha = 2/3;
	    beta  = -1/3;
	}else if (field=="Bnd"){
	    f.deltas = e.deltas^(-1) *sqrt(1-e.deltas); # e0^(-1)  * sqrt(2*(1-delta))
        alpha = 1;
        beta  = 1/2;
    }else{
        fprintf('Field %s not recognized.\n', field)
        stop("Invalid input!");

    }
	
	
	s = f.deltas * e.deltas;
	
    result <- list(eta=f.deltas,alpha=alpha,beta=beta, slope=s);
	
	return(result)
	
}
