function Phi = walsh( m )
%
N=2^m;
A=zeros(N,m);
count=0;
for c=N:2*N-1
    count=count+1;
    bin=dec2bin(c);
    bin=bin(2:m+1);
    
    A(count,:)=bin;
end
%Phi=(-1).^(A*A');
Phi=mod(A*A',2);
end

