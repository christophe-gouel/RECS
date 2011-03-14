function [X,grid] = spblkdiag(X,grid,output)
% SPBLKDIAG Computes a sparse block diagonal matrix
%
% Arguments:
%
%  X - pxqxn array, a list of n pxq arrays corresponding to a single block on
%  the diagonal.
%
%  grid (optional) - row and column indices for building the sparse matrix
%
%  output (optional, default: 1) - indicates if the matrix has to be calculated
%  or not (0 or 1)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin < 3
  output = 1;
  if nargin < 2
    grid = [];
  end
end

[p,q,n] = size(X);
% Generate the row and column indices for building the sparse matrix
if isempty(grid) || ~isequal(size(grid),[n*p*q 2])
  P = (1:p)';
  Q = (1:q)';
  grid_sparse = [P(P(:,ones(q,1)),:) Q(Q(:,ones(p,1))',:)];

  grid        = zeros(n*p*q,2);
  incr        = [p q];
  incr        = incr(ones(size(grid_sparse,1),1),:);
  for i=1:n
    grid(((i-1)*p*q+1):(i*p*q),:) = grid_sparse+(i-1)*incr;
  end
end

if output
  X = sparse(grid(:,1),...
             grid(:,2),...
             X(:),...
             grid(end,1),...
             grid(end,2));
else
  X = [];
end
