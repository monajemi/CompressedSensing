% This utility multiplies an nx1 vector x by an mxn binary matrix A
% of order q^2 x n , where q is any integer (though in compressed
% sensing applications q is a prime number). The output of this
% function is "y=Ax".
% The matrix A is represented by an equivalent matrix B of order qxn 
% which stores the indices of ones in each block of each column.

function y = Exp_mult(B,x) 

S = size(B) ;
q = S(1) ; n = S(2) ; m = q^2  ;

% Initialize
y = 0*ones(m,1) ;

for j=1:n
for l=1:m
b = floor((l-1)/q) + 1 ;
c = l - 1 - (b-1)*q ;
if B(b,j) == c
y(l) = y(l) + x(j) ;
end
end
end

end
