function bin_Phi = gray( Phi )
%GRAY
R=real(Phi);
S=imag(Phi);
bin_Phi=[R-S;R+S];
end

