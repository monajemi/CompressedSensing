function A = buildGaussian(n,N,field)
%
% Author: Hatef Monajemi
% Aug 1, 2018

if nargin < 3,
  field = 'R';
end

switch field,
   case {'O'}
   		A = buildOctonionGaussian(n,N,field);
   case {'Q'}
   		A = buildQuaternionGaussian(n,N,field);
   case {'C'}
        A = buildComplexGaussian(n,N,field);
   case {'R2'}
        A = buildGaussianPairs(n,N,field);
   case {'R','R+','Pos','Bnd'}
        A = buildGaussianPurelyReal(n,N,field);
end













function A = buildOctonionGaussian(n,N,field)

%% Reference Hatef's note:
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

A = zeros(8*n, 8*N);

%% draw each z from gaussian N(0, 1/8)

z1 = randn(n,N)/sqrt(8);
z2 = randn(n,N)/sqrt(8);
z3 = randn(n,N)/sqrt(8);
z4 = randn(n,N)/sqrt(8);

z5 = randn(n,N)/sqrt(8);
z6 = randn(n,N)/sqrt(8);
z7 = randn(n,N)/sqrt(8);
z8 = randn(n,N)/sqrt(8);


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
		 
	A(8*(i-1)+1:8*i, 8*(j-1)+1:8*j) = Z;
  end
end

end









function A = buildQuaternionGaussian(n,N,field)

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

A = zeros(4*n, 4*N);

z1 = randn(n,N)/sqrt(4);
z2 = randn(n,N)/sqrt(4);
z3 = randn(n,N)/sqrt(4);
z4 = randn(n,N)/sqrt(4);



% fill out the matrix will little 4x4 matices
for i = 1:n
  for j = 1:N
  	Z = [z1(i,j) -z2(i,j) -z3(i,j) -z4(i,j); ...
  		 z2(i,j)  z1(i,j) -z4(i,j)  z3(i,j); ...
  	     z3(i,j)  z4(i,j)  z1(i,j) -z2(i,j); ...
  	     z4(i,j) -z3(i,j)  z2(i,j)  z1(i,j)]; 
	A(4*(i-1)+1:4*i, 4*(j-1)+1:4*j) = Z;
  end
end


end












function A = buildGaussianPurelyReal(n,N,field)
A = randn(n,N);
end






function A = buildComplexGaussian(n,N,field)

A = randn(n,N)/sqrt(2) + sqrt(-1) * randn(n,N)/sqrt(2) ;

end



function At = buildGaussianPairs(n,N,field)
if nargin < 3,
  field ='R2';
end

A1 = randn(n,N)/sqrt(2) + sqrt(-1) * randn(n,N)/sqrt(2) ;

A = [  real(A1)  -imag(A1) ;
       imag(A1)   real(A1)];

At = zeros(size(A));
for j=1:N,
    At(:,2*j-1) = A(:,j  );
    At(:,2*j  ) = A(:,j+N);
end

end
end




  
