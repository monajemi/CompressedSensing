% building a DG frame explicitly
% Author: Hatef Monajemi
% Date: 2012 May 29 
% 
% flag: can be 'REAL' , 'IMAG'
% k : number of nnz
% m : 2^{m} for Cplex is the number of rows (fixed)
%     2^(m+1) REAL   
% delta : n/N (1/delta has to be an integer)
% r: DG(m,r) refer to Sina Jafarpour's thesis


function A = buildDG(m, r, delta, AA, flag)
%u: xPx^T (P: Kerdoc) b: xb^T (b: vector)
 ratio = int32(1/delta);
 % sanity check
	if (ratio > 2^ (m*(r+1)) )
	 error('ratio cannot be more than 2^{m(r+1)}')
	end
	
	
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
otherwise
warning('unknown flag: must be REAL or IMAG') 
end


	
	% build frame
	A = [];
	H = hadamard(n);
	for i = 1:ratio
		A = [A , diag(AA.u(:,i)) * H];
	end
	A = A ./ sqrt(n);  % normalization
	
end
