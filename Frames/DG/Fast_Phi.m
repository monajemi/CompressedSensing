function f = Fast_Phi( alpha,  u )
%u: xPx^T (P: Kerdoc) b: xb^T (b: vector)
N=size(u,1);
C=size(alpha,1);
ratio=C/N;
alpha_temp=reshape(alpha,C/ratio,ratio);
alpha_temp=fht(alpha_temp);
f=sum(u.*alpha_temp,2);
f=f/sqrt(N);
end

