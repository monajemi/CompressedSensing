function alpha = PhiT( f,  u,b ,range)
%u: xPx^T (P: Kerdoc) b: xb^T (b: vector)
N=size(b,2);
%%%u=u(2:N,:);
%%%b=b(2:N,:);
ratio=size(u,2);
supp=find(f);
s_supp=size(supp,1);
u=u(supp,:); b=b(supp,:);
supp_size=length(supp);
u=repmat(u',1,N);
u=reshape(u',s_supp,N*ratio)';
b=repmat(b',ratio,1);
if(nargin<4)
    alpha=sqrt(-1).^(2*b-u)*f(f~=0);
else
   alpha=sqrt(-1).^(2*b(range,:)-u(range,:))*f(f~=0); 
end
end


