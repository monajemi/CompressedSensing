function Phi = gen_frame( m,ratio,vals )
if(nargin<3)
    vals=0;
end
N=2^(m+1); 
C=N*ratio/2;
u=gsm(m,0);
u=u(:,1:ratio);
i=sqrt(-1);
u=i.^u;
u=gray(u);
b=walsh((m+1));
if(vals==1)
    vals_b=repmat(b,1,ratio);
    vals_u=repmat(u,N,1);
    vals_u=reshape(vals_u,N,ratio*N);
    Phi.vals=mod(vals_u+2*vals_b,4);
   % Phi.vals=Phi.vals(2:N,:);
    %maybe we do not need this line:
    Phi.vals=sqrt(-1).^Phi.vals;
end

%Phi.b=b;
Phi.u=u;
end

