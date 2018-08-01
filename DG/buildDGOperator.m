% Operatorized DG frame 
% Author: Hatef Monajemi, Sina Jafarpour
% Date: 2012 May 29 
% 
% flag: can be 'REAL' , 'IMAG'
% k : number of nnz
% m : 2^{m} for Cplex is the number of rows (fixed)
%     2^(m+1) REAL   
% delta : n/N (1/delta has to be an integer)
% r: DG(m,r) refer to Sina's thesis


function A = buildDGOperator(m, r, delta, AA, flag)


%check1 ==  
ratio = int32(1/delta);
% define a DG operator
switch flag
  
  case {'IMAG','C'} 
  	n = 2^m;
  	N = ratio * n;

  case {'REAL','R'}
    n = 2^(m+1);
    N = ratio * n;	
   
   case 'Pos'
	n = 2^(m+1);
	N = ratio * n;
   case 'Bnd'
	n = 2^(m+1);
	N = ratio * n;
end

	A = operator(@(w,mode)MatrixVectorMult(m,r, delta,AA,w,mode),n,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Private functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = MatrixVectorMult(m, r, delta, AA,w,mode)
% Spike+Sines operator.
% n : # of rows
% If mode == 1, returns  v =  [I F]  *w;
% if mode == 2, returns  v = [I F]' * w,

	if mode == 1
    ratio = int32(1/delta);
    Phi= @(x) Fast_Phi(x,AA.u);
    v=Phi(w);
    
	elseif mode == 2
   
    ratio = int32(1/delta);
    Phi = @(x) Fast_PhiT(x,AA.u);
    v=Phi(w);
    
    end
    
	end


end
