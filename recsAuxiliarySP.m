function [interp,xa,za] = recsAuxiliarySP(model,interp)
% RECSAUXILIARYSP calculates auxiliary variables not included in the core model
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z)
% For the definition of MODEL and INTERP structure, see recsSolveREE
% documentation. RECSAUXILIARY returns the interpolation structure INTERP, in
% which the coefficients matrices cxa (for auxiliary variables) and cza (for
% auxiliary expectations) have been added as new fields.
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z,XA) uses XA as first guess for the
% value of the auxiliary variables on the grid. This is only useful when
% auxiliary variables are function of auxiliary expectations. In this case, it
% can reduce a lot the time to reach the solution to start from a good first
% guess.
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z,XA,OPTIONS) solves for auxiliary
% variables with the parameters defined by the structure OPTIONS. The fields of
% the structure are
%    eqsolver         : 'fsolve', 'krylov' (default), 'lmmcp', 'ncpsolve',
%                       'path' or 'sa'
%    extrapolate      : 1 or 2 if extrapolation is allowed outside the
%                       interpolation space, 0 or -1 to forbid it (default: 1).
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [INTERP,XA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of the
% auxiliary variables on the grid.
%
% [INTERP,XA,ZA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of
% the auxiliary expectations on the grid.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2017 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
    
%% Initialization
functions = model.functions;
shocks    = model.shocks;
params    = model.params;

cX     = interp.cx;
fspace = interp.fspace;
Phi    = interp.Phi;
s      = interp.s;
X      = interp.X;

k = cell(nperiods,1); n = k; za = k; ss = k; ee = k; xx = k; snext = k;
Phinext = k; ind = k; z = k; cxa = k; xa = k;

atol      = sqrt(eps);
maxit     = 1000;
showiters = true;

dim       = cell2mat(model.dim{:});
p         = dim(:,3);
dima      = cell2mat(model.dima{:});
pa        = dima(:,2);
inext     = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);

for i=1:nperiods
  za{i}           = zeros(size(s{i},1),pa{i});
  k{i}            = size(shocks{i}.e,1); 
  n{i}            = size(s{i},1);
  ind{i}          = (1:n{i});
  ind{i}          = ind{i}(ones(1,k{i}),:);
  ss{i}           = s{i}(ind,:);
  ee{i}           = shocks{i}.e(repmat(1:k{i},1,n{i}),:);
  xx{i}           = X{i}(ind,:);
  snext{i}        = functions(i).g(ss{i},xx{i},ee{i},params,struct('F',1,'Js',0,'Jx',0));
  Phinext{i}      = funbas(fspace{inext(i)},snext{i});
  [LBnext,UBnext] = functions(inext(i)).b(snext{i},params);
  xnext           = min(max(Phinext{inext(i)}*cX{inext(i)},LBnext),UBnext);
  h               = functions(i).h(ss{i},xx{i},ee{i},snext{i},xnext{i},params);
  z{i}            = reshape(shocks{i}.w'*reshape(h,k{i},n{i}*p(i)),n{i},p(i));
end


%% Calculate the auxiliary variables
for i=1:nperiods
  xa{i} = functions(i).fa(s{i},X{i},z{i},za{i},params);
  cxa{i} = funfitxy(fspace{i},Phi{i},xa{i});
end
if any(pa>0)
  %% With auxiliary expectations function
  cnrm     = 1;
  it       = 0;
  vec      = @(x) cell2mat(cellfun(@(z) z(:),x,'UniformOutput',false));
  
  if showiters
    fprintf(1,'Successive approximation\n');
    fprintf(1,'   Iter\tResidual\n');
  end
  
  while(cnrm > atol && it < maxit)
    it    = it+1;
    cxaold = cxa;
    
    for i=nperiods:-1:1
      xanext = Phinext{inext(i)}*cxa{inext(i)};
      ha     = functions(i).ha(ss{i},xx{i},xa{i}(ind{i},:),ee{i},snext{i},xnext{i},xanext,params);
      za{i}  = reshape(shocks{i}.w'*reshape(ha,k{i},n{i}*pa(i)),n{i},pa(i));
      xa     = functions(i).fa(s{i},X{i},z{i},za{i},params);
      cxa{i} = funfitxy(fspace{i},Phi{i},xa{i});
    end % for i=nperiods:-1:1
 
  cnrm = norm(vec(cxa)-vec(cxaold));
  if showiters, fprintf(1,'%7i\t%8.2E\n',it,cnrm); end
  end % while(cnrm > atol && it < maxit)
end

%%
interp.cxa = cxa;
