# PTtools
These are codes for analytical predictions of compressed sensing Phase Transitions and various function for the analysis of PT datasets. Not everything is currently made available. I will make more functions available over time.

To install the R package of PT from github:
```
> install.packages("devtools")
> library(devtools)
> install_github("monajemi/CompressedSensing", subdir="PTtools/R/PT")
> library(PT)

> predictPT(1/2,'R')
[1] 0.1928448
```

There is also a MATLAB version of this prediction.

```
>> predictPT(1/2,'R')
```


# Frames
These are codes for construction of various frames that I have used in my publications 

## Main function
```
function  A = buildFrame(n, N, ensembleCode,fieldCode,qOperator,varargin);
% Function: buildFrame.m
% builds desired frame according to information given as input.
%
% input:
% n: # rows
% N: # columns
% ensembleCode: the code for the desired ensemble e.g., 'DG'
% fieldCode: 'R','C','R+'
% qOperator: if set to 1, then an operator A is built for fast calculation of Ax and A'x.
%            if set to 0 (default), frame is explicitly built
% varargin:
% LDPC and devore/array_pc: 
%   varargin{1}:column degree
%   varargin{2}:sglConcent (0: strict regularity, 1:best-effort)
```

## Examples: 
```
% build a real USE matrix
A = buildFrame(16,64,'USE','R');

% build a Complex DG frame of size 32x128
A = buildFrame(32,32*4,'DG','C');


% build 'non-strictly' regular (best effort) LDPC matrices of degree 3
A = buildFrame(504,1008, 'LDPC','R',[],3,1);
A = buildFrame(504,1008, 'LDPC','R',[],3);


% build 'strictly' regular LDPC matrices of degree 3
A = buildFrame(504,1008, 'LDPC','R',[],3,0);

% build a Devore matrix of size (q^2)x20 with coulmn degree q=3 
A = buildFrame(9,20,'devore',[],[],3)

% build an array parity check matrix of size (j*q)xq^2 for j=3 and q=5;
A = buildFrame(5*3,5^2,'array-pc',[],[],3)

```


## References: 

"Deterministic matrices matching the compressed sensing phase transitions of Gaussian random matrices", Monajemi et al. 2013. [Download](http://www.pnas.org/content/110/4/1181)

"Compressed Sensing Using Binary Matrices of Nearly Optimal Dimensions", Mahsa Lotfi and Mathukumalli Vidyasagar, [arXiv: 1808.03001](https://arxiv.org/abs/1808.03001)
