%%This function recieves the matrix of indices, "B", as the input and gives us
%%the measurement map, "A" which is the binary version of the DeVore's
%%construction. "A" is (q^2)xn.
function A = B_to_A(B)

S = size(B) ;
q = S(1) ;
n = S(2) ;
m = q^2 ;
A = 0*ones(m,n) ;
for j=1:n
for l=1:m
b = floor((l-1)/q) + 1 ;
c = l - (b-1)*q ;
if B(b,j) == c-1
A(l,j) = 1 ;
end
end
end

end
