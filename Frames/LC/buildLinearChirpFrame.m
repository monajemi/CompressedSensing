function A = buildLinearChirpFrame(n,N,field)
%  A = buildLinearChirpFrame(n,N,field)
% Inputs
%   n     number of samples
%   N     number of columns, N = n*L  L integral number of  sub-bases 1 < L <=n
%   field 'Bnd','Pos','R','C'
%
% Outputs
%   A     complex(n,N)
%         N == n*L
%
% Authors: David Donoho, Hatef Monajemi

if nargin < 3,
   field = 'R';
end

if ~isprime(n),
   error('PAST YOUR PRIME','Sorry, n=%i is not prime',n);
end

if abs(n*(N/n - round(N/n))) < .5,
  L = round(N/n);
else
  L = ceil(N/n);
end


A0 = buildComplexLinearChirpFrame(n,L);
switch field,
 case {'R','R+', 'Pos','Bnd'},
  A = embedCplexAsRealInstance(A0);
  N0  = size(A0,2);
  N1  = size(A,2);
  ex = N1/N0;
  if N*ex<N1,
     A = trimMatrix(A,n*ex,N*ex,N1,field);
  end
 case  'R2',
    A = embedCplexAsRealInstance(A0);
 case  'C',
    A = A0;
end

end

function A = buildComplexLinearChirpFrame(n,L)
%  A = buildComplexLinearChirpFrame(n,L)
% Inputs
%   n     number of samples
%   L  number of sub-bases 1 < L <=n
% Outputs
%   A     complex(n,N)
%         N == n*L
%
    if nargin < 2,
       L = n;
    end
    %
    t = (0:(n-1))';
    step = floor(n/L);
    N = L*n;
    A = zeros(n,N);
    F = fft(eye(n)) ./sqrt(n);
    for ell=1:L,
        g = exp(sqrt(-1).*2*pi/n.* (ell-1) .*t.^2);
        D = diag(g);
        this_inx = ((ell-1)*n+1):(ell*n);
        A(:,this_inx) = D * F;
    end

end

function A0 = trimMatrix(A,n,N,N1,field)
  A0 = A;
  if(N < N1),
     [~,j0] = sort(rand(N1,1));
     switch field,
      case {'R','R+','Pos','Bnd','C'}
          A0 = A(:,j0(1:N));
      case {'R2'}
          sel = FALSE(2*N1,1);
          for i=1:N,
             k0 = 2*j0(i)-1;
             sel(k0)  = TRUE;
             sel(k0+1) = TRUE;
          end
          A0 = A(:,sel);
     end
  else
    error('Ncredible','n(=%i),N(=%i),N1(=%i)',n,N,N1)
  end
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







%end
