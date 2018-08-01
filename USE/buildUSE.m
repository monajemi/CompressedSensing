function A = buildUSE(n,N,field)
if nargin < 3,
  field = 'R';
end

switch field,
   case {'O'}
   		A = buildOctonionUSE(n,N,field);  
   case {'Q'}
   		A = buildQuaternionUSE(n,N,field);  
   case {'C'}
        A = buildComplexUSE(n,N,field);
   case {'R','R+','Pos','Bnd'}
        A = buildUSEPurelyReal(n,N,field);
end













function A = buildOctonionUSE(n,N,field)

%% Reference Hatef Monajemi's thesis:
%% An Octonion number is represented by a 8x8 matrix
%% 
%% z = z
%% 
%% Z = See below for matrix representaion 
%% This representation is not commutative. 
%%
%% z . q = [Z] q (simple matrix vector product on the right)
%%
%%

% allocate memory for unnormalized sensing matrix A0;

A0 = zeros(8*n, 8*N);

z1 = randn(n,N);
z2 = randn(n,N);
z3 = randn(n,N);
z4 = randn(n,N);

z5 = randn(n,N);
z6 = randn(n,N);
z7 = randn(n,N);
z8 = randn(n,N);


% fill out the matrix will little 8x8 matices
% check http://arxiv.org/pdf/math/0003166v2.pdf
for i = 1:n
  for j = 1:N
  
  	Z11 =[z1(i,j) -z2(i,j) -z3(i,j) -z4(i,j) ; ...
  		 z2(i,j)  z1(i,j) -z4(i,j)  z3(i,j) ; ...
  	     z3(i,j)  z4(i,j)  z1(i,j) -z2(i,j) ; ...
  	     z4(i,j) -z3(i,j)  z2(i,j)  z1(i,j)];
  	
  	Z12 = [-z5(i,j) -z6(i,j)  -z7(i,j)  -z8(i,j); ...
  		 -z6(i,j)  z5(i,j)   z8(i,j)  -z7(i,j); ...
  	     -z7(i,j) -z8(i,j)   z5(i,j)   z6(i,j); ...
  	     -z8(i,j)  z7(i,j)  -z6(i,j)   z5(i,j)]; 

	Z21 = [z5(i,j) z6(i,j)  z7(i,j)  z8(i,j); ...
  		   z6(i,j)  -z5(i,j)   z8(i,j)  -z7(i,j); ...
  	       z7(i,j) -z8(i,j)   -z5(i,j)   z6(i,j); ...
  	       z8(i,j)  z7(i,j)  -z6(i,j)   -z5(i,j)]; 

	Z22 =[z1(i,j) -z2(i,j) -z3(i,j) -z4(i,j) ; ...
  		  z2(i,j)   z1(i,j)   z4(i,j) -z3(i,j) ; ...
  	      z3(i,j)  -z4(i,j)  z1(i,j)  z2(i,j) ; ...
  	      z4(i,j)   z3(i,j) -z2(i,j)  z1(i,j)];
  	

	Z = [Z11, Z12;...
		 Z21, Z22];
		 
	A0(8*(i-1)+1:8*i, 8*(j-1)+1:8*j) = Z;
  end
end

D0 = diag(A0'*A0);
A  = A0 * diag(1./sqrt(D0));

end









function A = buildQuaternionUSE(n,N,field)

%% Reference Hatef's note:
%% A quaternion number is represented by a 4x4 matrix
%% 
%% z = z1+ z2 i + z3 j + z4 k 
%% 
%% Z = [z1 -z2 -z3 -z4; z2, z1, -z4, z3; z3, z4, z1, -z2; z4, -z3, z2, z1] 
%% This representation is not commutative. This corresponds to Ell's approach.
%%
%% z . q = [Z] q (simple matrix vector product on the right)
%%
%%

% allocate memory for unnormalized sensing matrix A0;

A0 = zeros(4*n, 4*N);

z1 = randn(n,N);
z2 = randn(n,N);
z3 = randn(n,N);
z4 = randn(n,N);



% fill out the matrix will little 4x4 matices
for i = 1:n
  for j = 1:N
  	Z = [z1(i,j) -z2(i,j) -z3(i,j) -z4(i,j); ...
  		 z2(i,j)  z1(i,j) -z4(i,j)  z3(i,j); ...
  	     z3(i,j)  z4(i,j)  z1(i,j) -z2(i,j); ...
  	     z4(i,j) -z3(i,j)  z2(i,j)  z1(i,j)]; 
	A0(4*(i-1)+1:4*i, 4*(j-1)+1:4*j) = Z;
  end
end

D0 = diag(A0'*A0);
A = A0 * diag(1./sqrt(D0));


end












function A = buildUSEPurelyReal(n,N,field)
A = randn(n,N);
for j=1:N,
  A(:,j) = A(:,j)./norm(A(:,j),'fro');
end

end






function A= buildComplexUSE(n,N,field)

A0 = randn(n,N) + sqrt(-1) * randn(n,N) ;
D0 = diag(A0'*A0);
A = A0 * diag(1./sqrt(D0));


end






end




  
