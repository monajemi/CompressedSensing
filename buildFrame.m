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
%
%
% Copyright 2012: Hatef Monajemi(monajemi@stanford.edu), David Donoho (donoho@stanford.edu)

% add path to all subdirectories
addpath(genpath('.'))


if nargin < 4, warning('the number of arguments must be at least 4'); end
if nargin < 5 || isempty(qOperator), qOperator = 0; end


       switch ensembleCode,
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
