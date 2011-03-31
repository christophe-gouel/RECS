function recsAccuracy(model,interp,s)
% RECSACCURACY Evaluates approximation accuracy
%
% RECSACCURACY(MODEL,INTERP,S)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

params = model.params;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func   = model.func;
else
  error('model.func must be either a string or a function handle')
end

cx     = interp.cx;
cz     = interp.cz;
fspace = interp.fspace;

disp('Accuracy of the solution');

[n,d,t] = size(s);
se      = permute(s,[1 3 2]);
se      = reshape(se,n*t,d);

[LB,UB] = func('b',se,[],[],[],[],[],params);
xe      = min(max(funeval(cx,fspace,se),LB),UB);
ze      = funeval(cz,fspace,se);
EE      = func('e',se,xe,ze,[],[],[],params);
lEE     = [log10(max(abs(EE)));
           log10(sum(abs(EE))/size(EE,1))];

disp(' Euler equation error (in log10)');
disp('    Max       Mean');
disp(lEE');

fe      = func('f',se,xe,ze,[],[],[],params);

Ef = min(max(-fe,LB-xe),UB-xe);
Ef = [max(abs(Ef));
      mean(abs(Ef))];

disp(' Equilibrium equation error (minmax formulation)');
disp('    Max       Mean');
disp(Ef');


