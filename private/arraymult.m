function c = arraymult(a,b,n,p,q,r)
% ARRAYMULT Computes matrix multiplication for 3-D arrays
%
% For a faster execution, it is better to install the compiled mex files from
% the CompEcon toolbox.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
c = zeros(n,p,r);
if r > 1
  for i = 1:p
    for k = 1:q
      c(:,i,:) = c(:,i,:) + a(:,i,k*ones(r,1)).*b(:,k,:);
    end
  end
else % r==1
  for k=1:q
    c = c + a(:,:,k).*b(:,k*ones(p,1));
  end
end
