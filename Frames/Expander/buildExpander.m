% the coo implementation of expander matrix

% Notes by Jared Tanner:
% A somewhat better implementation would be to use sparse matrix format such as coo,
% and also not to regenerate the random permutation for each column.  It is preferable to
% generate the matrix as in the attached, which generates a random permutation of m, uses the
% first d of these as the row index for the first column, uses the next d of them as the row index
% for the second column, etc.. until needing to generate a new random permutation of m for more
% rows.  This ensures that the first floor(m/d) columns are disjoint which isn’t true for the random
% construction
% Author: Hatef Monajemi, Jared Tanner
% Date: Sept 22,2018

function A = buildExpander(n,N,field,d)
if nargin < 3,
field = 'R';
d = 7;
end

% nnz is number of nonzeros
% nperm is number of permutations
% needed

nnz = d*N;

col_per_perm=floor(n/d);
nperms=ceil(N/col_per_perm);

val=ones(nnz,1);
col=zeros(nnz,1);
row=zeros(nnz,1);

col_ind=1;
for j=1:nperms
for k=1:col_per_perm
if col_ind<=N
col(1+(col_ind-1)*d:col_ind*d)=col_ind*ones(d,1);
end
col_ind=col_ind+1;
end
end

for j=1:nperms-1
row(1+(j-1)*n:j*n)=randperm(n);
end
row(j*n+1:end)=randperm(length(row(j*n+1:end)));

A = sparse(row,col,val);


end


