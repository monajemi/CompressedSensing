% This function generates a matrix B of dimensions q by n where q is a prime,
% and each element of B belongs to the prime number field Z/(q).

function B = DeVore_Exp(q , n)

if isprime(q)==0
display('Error: input must be a prime number')
return
end

if n >= q^3+1
display('Error: n must be no larger than q^3')
return
end

% Initialize B matrix to all zeros.

B = zeros(q,n) ;

% Now we construct all possible polynomials of degree 2 and evaluate them at all points in Z/(q).

count = 1 ;
for i=0:q-1
    for j=0:q-1
        for k=0:q-1
			coef = [ i j k ] ;
			L = 0:1:q-1 ;
			B(:,count) = mod(polyval(coef,L) , q);
			count = count+1 ;
			if count == n+1
				break ;
			end
        end
		if count == n+1
			break ;
		end
	end
	if count == n+1
		break ;
	end
end

end


