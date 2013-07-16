function Diagnostics(X,F,J)
% DIAGNOSTICS Provides diagnostics information about a problem of nonlinear equations

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[xmax,ix] = max(abs(X));
[fmax,iF] = max(abs(F));

iJ           = NaN(2,1);
[Jmax,tmp]   = max(abs(J));
[Jmax,iJ(2)] = max(Jmax);
Jmax         = full(Jmax);
iJ(1)        = tmp(iJ(2));

Jstat  = NaN(4,1);
iJstat = zeros(4,1);
for i=1:size(J,1)
  tmp = norm(J(i,:),2);
  Jstat(1) = max(Jstat(1),tmp);
  if Jstat(1)==tmp, iJstat(1) = i; end
  Jstat(2) = min(Jstat(2),tmp);
  if Jstat(2)==tmp, iJstat(2) = i; end
  tmp = norm(J(:,i),2);
  Jstat(3) = max(Jstat(3),tmp);
  if Jstat(3)==tmp, iJstat(3) = i; end
  Jstat(4) = min(Jstat(4),tmp);
  if Jstat(4)==tmp, iJstat(4) = i; end
end

disp(repmat('_',1,60))
fprintf(1,'Diagnostic Information\n\n')

fprintf(1,'Number of variables: %d\n\n',length(X))

fprintf(1,'%d nonzeros Jacobian elements, %3.2f%% dense \n\n',nnz(J)*[1 100/numel(J)])

disp('Point statistics')
fprintf(1,'Maximum of X. . . . . . . . . .  %6.4e var: (x[%5d])\n',[xmax ix])
fprintf(1,'Maximum of F. . . . . . . . . .  %6.4e eqn: (f[%5d])\n',[fmax iF])
fprintf(1,'Maximum of Grad F . . . . . . .  %6.4e eqn: (f[%5d])\n',[Jmax iJ(1)])
fprintf(1,'                                            var: (x[%5d])\n',iJ(2))

disp('Jacobian norm statistics')
fprintf(1,'Maximum Row Norm. . . . . . . .  %6.4e eqn: (f[%5d])\n',[Jstat(1) iJstat(1)])
fprintf(1,'Minimum Row Norm. . . . . . . .  %6.4e eqn: (f[%5d])\n',[Jstat(2) iJstat(2)])
fprintf(1,'Maximum Column Norm . . . . . .  %6.4e eqn: (x[%5d])\n',[Jstat(3) iJstat(3)])
fprintf(1,'Minimum Column Norm . . . . . .  %6.4e eqn: (x[%5d])\n',[Jstat(4) iJstat(4)])

disp(repmat('_',1,60))