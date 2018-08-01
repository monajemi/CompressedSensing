function f = Phi( alpha,  u,b )
%u: xPx^T (P: Kerdoc) b: xb^T (b: vector)
N=size(b,2);
supp=find(alpha)-1;
f=sqrt(-1).^(2*b(:,1+mod(supp,N))+u(:,1+floor(supp/N)))*alpha(alpha~=0);
%%%f=f(2:N,:);
end

