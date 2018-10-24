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
% LDPC: 
%   varargin{1}:column degree
%   varargin{2}:sglConcent (0: strict regularity, 1:best-effort)
%
%
% Copyright 2012: Hatef Monajemi(monajemi@stanford.edu), David Donoho (donoho@stanford.edu)
% Relevant Citation: "Deterministic matrices matching the compressed sensing phase transitions of Gaussian random matrices.", 
% Monajemi et al., 2013 http://www.pnas.org/content/110/4/1181 
%
%
% Examples: 
%% build a real USE matrix
% A = buildFrame(16,64,'USE','R');
%
%% build Exapander of col degree 7
% A = buildFrame(100,300,'Expander','R',[],7);
%
%% build 'non-strictly' regular (best effort) LDPC matrices of degree 3
% A = buildFrame(504,1008, 'LDPC','R',[],3,1);
% A = buildFrame(504,1008, 'LDPC','R',[],3);
%
%% build 'strictly' regular LDPC matrices of degree 3
% A = buildFrame(504,1008, 'LDPC','R',[],3,0);



% add path to all subdirectories
addpath(genpath('.'))


if nargin < 4, warning('the number of arguments must be at least 4'); end
if nargin < 5 || isempty(qOperator), qOperator = 0; end


       switch ensembleCode,

                case  'array-pc'
                    if length(varargin) < 1
                        error('please set the column degree of the array parity check matrix to more than 0 and less than sqrt(n)')
                    else
                        A = Array_Parity_Check_Matrix(N,varargin{1});
                    end
                
                case  'devore'
                    if isprime(varargin)==0 %%column weight must be a prime number, n=(coulmn_weight)^2 and N<=(coulmn_weight)^3 
                        error('please set the column degree of the devore matrix to a prime number that is more than 1')
                    else
                         if varargin>=sqrt(N) or varargin< nthroot(N,3)
                              error('please set the column degree of the devore matrix to less than square root of N and greater than or equal to cube root of N')
                         else
                              B = DeVore_Exp(varargin{1},N);
                              A = B_to_A(B);
                         end
                    end
                
                
                case  'Gaussian'
                    if qOperator== 1
                        error('Operator not supported for this case')
                    else
                        A = buildGaussian(n,N,fieldCode);
                    end


                case 'Expander'
                if length(varargin) < 1
                    error('please enter the column degree of the Exapnder matrix')
                end
                A = buildExpander(n,N,fieldCode,varargin{1});


                case 'LDPC'
                    if length(varargin) < 1
                        error('please enter the column degree of the LDPC matrix')
                    end
                    
                    
                    % LDPC needs mex on cluster
                    mex LDPC/GenerateLDPC_mex.C LDPC/BigGirth.C LDPC/CyclesOfGraph.C LDPC/Random.C LDPC/GenerateLDPC.C

                    

                    if qOperator== 1
                        error('Not Available');
                    else
                        % Only real matrix
                        if (length(varargin) == 2)
                            if ~ ( varargin{2} == 0 || varargin{2} == 1)
                                error('sglConcent must be 0 or 1. 0: strictly regular, 1:best effort')
                            end
                            A = buildLDPCmatrix(n,N,varargin{1},varargin{2});
                        else
                            A = buildLDPCmatrix(n,N,varargin{1});
                        end    
                    end

                case  'LC'
                    if qOperator== 1
                    error('Not yet implemented');
                    else
                    A = buildLinearChirpFrame(n,N,fieldCode);
                    end


                case  'GF'
                    if qOperator== 1
                    error('Not yet implemented');
                    else
                        A = buildGrassmannianFrame(n,N,fieldCode);
                    end

                case  'PETF'
                    if qOperator== 1
                        error('Operator not supported for this case');
                    else
                        A  = buildPaleyETF(n,N,fieldCode);
                    end
                case  'USE'
                    if qOperator== 1
                       error('Operator not supported for this case')
                    else
                       A = buildUSE(n,N,fieldCode);
                    end  

                
                
                case  'DG'

                mex DG/gsm_mex.cpp

                % determin mDG, and AA
                if( strcmp(fieldCode, 'C') )
                    % complex
                    mDG = log2(n);
                    if (mod(mDG,2) ~= 1)
                    error('DG(m,r) is valid for odd m');
                    end
                    AA = gen_frame_cplex(mDG, N/n, 0);

                else 
                    % real
                    mDG = log2(n)-1;
                    if (mod(mDG,2) ~= 1)
                    error('DG(m,r) is valid for odd m');
                    end
                    AA = gen_frame(mDG, N/n, 0);
                end


                if qOperator== 1
                    A = buildDGOperator(mDG, 0, n/N, AA,fieldCode);
                else
                    A = buildDG(mDG, 0, n/N, AA,fieldCode);
                end

                otherwise
                    error('ensemble %s is not recognized.', ensembleCode)
      end
end
