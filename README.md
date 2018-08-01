# PTpredictions
These are codes for analytical predictions of Phase Transition in Compressed Sensing

To install the R package of PT from github:
```
> install.packages("devtools")
> library(devtools)
> install_github("monajemi/CompressedSensing", subdir="PTtools/R/PT")
> library(PT)

> predictPT(1/2,'R')
[1] 0.1928448
```

# Frames
These are codes for construction of various frames that I have used in my publications 

# Example: 
```
% build a real USE matrix
A = buildFrame(16,64,'USE','R');

% build a Complex DG frame of size 32x1288
A = buildFrame(32,32*4,'DG','C');
```


## References: 

"Deterministic matrices matching the compressed sensing phase transitions of Gaussian random matrices", Monajemi et al. 2013
