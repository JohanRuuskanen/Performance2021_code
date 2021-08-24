function S=krons(A,B)
% S=KRONS(A,B)
% Kronecker sum of matrices A and B
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
S=kron(A,eye(size(B)))+kron(eye(size(A)),B);
end