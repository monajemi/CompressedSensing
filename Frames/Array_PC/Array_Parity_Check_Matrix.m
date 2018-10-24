%%This file contains implementation of the Array code parity check matrix
%%which has the girth equal to 6.The matrix has the dimesionality equal to 
%%mx(q^2) in which m=jq. Parameter "j" shows the column-weight of the matrix.


function H = Array_Parity_Check_Matrix(n,j) 
q=sqrt(n);    
I=eye(q);
P=zeros(q,q);
m=j*q;
H=zeros(m,n);
C=cell(j-1,q-1);
%%%find the permutation matrix P for I
for i=1:q
  [i1,i2]=find(I(i,:)==1);%%% i2 determines the index of 1 in each row i of I
  if mod(i2+(q-1),q)~=0
  i3=mod(i2+(q-1),q);
  P(i,i3)=1;  
  else
  P(i,q)=1; 
  end
end
%%%Now let's construct the H matrix
parfor k1=1:(j-1)%%%each row
    for k2=1:(q-1)
        
        C(k1,k2)={P^(k1*k2)};
    end
end
Partial=cell2mat(C);%%%this gives the part of the H that has the permutations in it
for k3=0:q-1
   H(1:q,k3*q+1:(k3+1)*q)=I;  
end
for k4=1:j-1
   H(k4*q+1:(k4+1)*q,1:q)=I;
end
H(q+1:j*q,q+1:q^2)=Partial;
end