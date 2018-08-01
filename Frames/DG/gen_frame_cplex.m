function Phi = gen_frame_cplex(m, ratio, vals )
if(nargin<3)
    vals=0;
end
N=2^(m); 
C=N*ratio;
u = gsm(m,0);
u=u(:,1:ratio);
i=sqrt(-1);
u=i.^u;

Phi.u = u;
end
