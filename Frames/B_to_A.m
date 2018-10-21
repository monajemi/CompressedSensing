function A = B_to_A(B)

S = size(B) ;
q = S(1) ;
n = S(2) ;
m = q^2 ;
A = zeros(m,n) ;
for j=1:n
	for l=1:m
		b = floor((l-1)/q) + 1 ;
		c = l - (b-1)*q +1;
		if B(b,j) == c
			A(l,j) = 1 ;
		end
	end
end

end
		