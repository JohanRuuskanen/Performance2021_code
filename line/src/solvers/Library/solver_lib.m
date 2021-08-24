function [QN,UN,RN,TN,CN,XN] = solver_lib(qn, options)
% [Q,U,R,T,C,X] = SOLVER_LIB(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%% generate local state spaces
M = qn.nstations;
K = qn.nclasses;
N = qn.njobs';
rt = qn.rt;
V = qn.visits;

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

line_error(mfilename,'This model or method is not supported yet. Returning with no result.');
end
