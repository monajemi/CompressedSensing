function A = buildGrassmannianFrame(n,N,field)
%  [y,A,x0] = buildGrassmannianFrame(k,n,N,field)
% Inputs
%   n     number of samples
%   N     number of columns, N = n*L  L integral number of sub-bases 1 < L <=n
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

L = round(N/n);
A0 = buildComplexGrassmannianFrame(n,L);

switch field,
 case {'R', 'R+', 'Pos','Bnd'},
  A = embedCplexAsRealInstance(A0);
  N0  = size(A0,2);
  N1  = size(A,2);
  ex = N1/N0;
  if N<L*n,
     A = trimMatrix(A,n*ex,N*ex,field);
  end
 case  'R2',
    A = embedCplexAsRealInstance(A0);
 case  'C',
    A = A0;
end

end

function A = buildComplexGrassmannianFrame(n,L)
%  A = buildComplexGrassmannianFrame(n,L)
% Inputs
%   n     number of samples
%   L     number of sub-bases 1 < L <=n
% Outputs
%   A     complex(n,N)
%         N == n*L
    if nargin < 2,
       L = n;
    end

    t = (0:(n-1))';
    alltop = exp( sqrt(-1) .* 2*pi * t.^3./n);
    step = floor(n/L);
    N = L*n;
    A = zeros(n,N);
    F = fft(eye(n)) ./sqrt(n);
    for ell=1:L,
        g = circshift(alltop,ell.*step);
        D = diag(g);
        this_inx = ((ell-1)*n+1):(ell*n);
        A(:,this_inx) = D * F;
    end

end

function x1 = circshift(x0,h)
    M = length(x0);
    if h > 0,
       x1 =  [ x0((M+1-h):M); x0(1:(M-h)) ];
    elseif h < 0,
       x1 =  [  x0((-h+1):M) ; x0(1:-h) ];
    else
       x1 = x0(:);
    end
end

function A0 = trimMatrix(A,n,N,field)
  A0 = A;
if(N < 2*n),
     [~,j0] = sort(rand(2*n,1));
     switch field,
      case {'R','R+','Pos','Bnd','C'}
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
       case {'R','R+','Pos','Bnd','C'}
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
