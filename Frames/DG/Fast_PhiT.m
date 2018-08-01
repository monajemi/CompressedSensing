function alpha = Fast_PhiT( f,  u,b )
%u: xPx^T (P: Kerdoc) b: xb^T (b: vector)


N=size(u,1);
ratio=size(u,2);
temp_f=repmat(f,1,ratio);
temp_f=temp_f.*conj(u);
temp_f=fht(temp_f);
alpha=reshape(temp_f,ratio*N,1);
alpha=(alpha/sqrt(N));
end


