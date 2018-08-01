function  spectrum= chirpogram(f,u,b ,m)
if(nargin<4)
    m=log2(size(b,1));
end
%f is N-1 times 1 measurement vector
d=deChirp(f,u);
%d is N times ratio chirpogram matrix [i^[-xPx^T].*f
spectrum=fht(d,m);
spectrum=reshape(abs(spectrum),size(u,2)*size(b,2) ,1);

end

