function recsAccuracy(model,interp,s)
% RECSACCURACY Evaluates approximation accuracy

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

func   = model.func;
params = model.params;

cx     = interp.cx;
cz     = interp.cz;
fspace = interp.fspace;

disp('Accuracy of the solution');

[n,d,t] = size(s);
se      = permute(s,[1 3 2]);
se      = reshape(se,n*t,d);

[LB,UB] = feval(func,'b',se,[],[],[],[],[],params);
xe      = min(max(funeval(cx,fspace,se),LB),UB);
ze      = funeval(cz,fspace,se);
EE      = feval(func,'e',se,xe,ze,[],[],[],params);
lEE     = [log10(max(abs(EE)));
           log10(sum(abs(EE))/size(EE,1))];

disp(' Euler equation error (in log10)');
disp('    Max       Mean');
disp(lEE');

fe      = feval(func,'f',se,xe,ze,[],[],[],params);

Ef = min(max(-fe,LB-xe),UB-xe);
Ef = [max(abs(Ef));
      mean(abs(Ef))];

disp(' Equilibrium equation error (minmax formulation)');
disp('    Max       Mean');
disp(Ef');


