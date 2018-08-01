function A = buildPaleyETF(n,N,field)
%
% Authors: David Donoho, Hatef Monajemi
%

if nargin < 3,
   field = 'R';
end


if (N < 2*n),
   p = 2*n-1;
elseif (n <= N/2),
   p = N-1;
end
%
if ~isprime(p),
    error('Prime Broker','PaleyETF needs p prime')
end

%  Paley TightFrame
switch field,
 case {'R','R+', 'Pos','Bnd'},
  A = realPaleyETF(p) ;
  A = trimMatrix(A,n,N,field);
 case {'C', 'R2'}
  A0 = PaleyETF(p);
  if strcmp(field,'R2'),
    A = embedCplexAsRealInstance(A0);
  else
    A = A0;
  end
end




function A0 = trimMatrix(A,n,N,field)
  A0 = A;
if(N < 2*n),
     [~,j0] = sort(rand(2*n,1));
     switch field,
      case {'R','R+','Pos','Bnd'}
          A0 = A(:,j0(1:N));
      case {'R2'}
          sel = FALSE(4*n,1);
          for i=1:N,
             k0 = 2*j0(i)-1;
             sel(k0)  = TRUE;
             sel(k0+1) = TRUE;
          end
          A0 = A(:,sel);
     end
elseif (2*n < N),
      switch field,
       case {'R','R+','Pos','Bnd'}
          [~,j0] = sort(rand(N/2,1));
          A0 = A(j0(1:n),:);
      case {'R2'} 
          sel = FALSE(N/2,1);     
          for i=1:n,
            k0 = 2*j0(i)-1;
            sel(k0)  = TRUE;
            sel(k0+1) = TRUE;
          end
          A0 = A(sel,:);
     end
  end
end

function rETF = realPaleyETF(p)
cETF = PaleyETF(p);
Gram = real(cETF'*cETF);
[V,D] = eig(Gram);
evalues = diag(D);
sel = real(evalues) > 1;
rV  = real(V(:,sel));
rETF = (rV * diag(sqrt(real(evalues(sel)))))';
end


function ETF = PaleyETF(p)
%
if ~isprime(p) || (rem(p,4) ~=1),
    error('PaleyETF needs p prime and == 1 mod 4')
end
M = (p-1)/2;              % number of quadratic residues
            
resi = rem((1:(p-1)).^2,p);
ures = unique(resi);
F = fft(eye(p));
H = F([1 1+ures],:);
D = diag([ sqrt(1/p) ; sqrt(2/p).*ones(M,1)]);
ETF = [D*H  [ 1 ; zeros(M,1)]]; 

end

function At = embedCplexAsRealInstance(A0)
%
[n,N] = size(A0);
A = [ real(A0)  -imag(A0) ;
      imag(A0)   real(A0)];

At = zeros(size(A));
for j=1:N,
    At(:,2*j-1) = A(:,j  );
    At(:,2*j  ) = A(:,j+N);
end

end

end


  
