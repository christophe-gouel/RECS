function [interp,xa,za,xass] = recsAuxiliarySP(model,interp)
% RECSAUXILIARYSP calculates auxiliary variables not included in the core model
%
% INTERP = RECSAUXILIARYSP(MODEL,INTERP)
% For the definition of MODEL and INTERP structure, see recsSolveREESP
% documentation. RECSAUXILIARYSP returns the interpolation structure INTERP, in
% which the coefficients matrices cxa (for auxiliary variables) and cza (for
% auxiliary expectations) have been added as new fields.
%
% [INTERP,XA] = RECSAUXILIARYSP(MODEL,INTERP) returns the value of the
% auxiliary variables on the grid.
%
% [INTERP,XA,ZA] = RECSAUXILIARYSP(MODEL,INTERP) returns the value of the
% auxiliary expectations on the grid.
%
% [INTERP,XA,ZA,XASS] = RECSAUXILIARYSP(MODEL,INTERP) returns the value of the
% auxiliary variables at the deterministic steady state.
%
% See also RECSSIMULSP, RECSSOLVEREESP.

% Copyright (C) 2011-2017 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
    
%% Initialization
functions = model.functions;
shocks    = model.shocks;
params    = model.params;
nperiods  = model.nperiods;

cX     = interp.cX;
fspace = interp.fspace;
Phi    = interp.Phi;
s      = interp.s;
X      = interp.X;

k = cell(nperiods,1); za = k; ss = k; ee = k; xx = k; snext = k;
Phinext = k; ind = k; z = k; cxa = k; xa = k; xnext = k; xass = k;
n = zeros(nperiods,1);

atol      = sqrt(eps);
atol      = 1E-9;
maxit     = 5E2;
showiters = 50;

p         = cell2mat(model.dim(:,3));
ma        = cell2mat(model.dima(:,1));
pa        = cell2mat(model.dima(:,2));
inext     = @(iperiod) (iperiod+1)*(iperiod<nperiods) + 1*(iperiod==nperiods);

for i=1:nperiods
  n(i)            = size(s{i},1);
  za{i}           = zeros(n(i),pa(i));
  k{i}            = size(shocks{i}.e,1); 
  ind{i}          = (1:n(i));
  ind{i}          = ind{i}(ones(1,k{i}),:);
  ss{i}           = s{i}(ind{i},:);
  ee{i}           = shocks{i}.e(repmat(1:k{i},1,n(i)),:);
  xx{i}           = X{i}(ind{i},:);
  snext{i}        = functions(i).g(ss{i},xx{i},ee{i},params,struct('F',1,'Js',0,'Jx',0));
  Phinext{i}      = funbas(fspace{inext(i)},snext{i});
  [LBnext,UBnext] = functions(inext(i)).b(snext{i},params);
  xnext{i}        = min(max(Phinext{i}*cX{inext(i)},LBnext),UBnext);
  h               = functions(i).h(ss{i},xx{i},ee{i},snext{i},xnext{i},params);
  z{i}            = reshape(shocks{i}.w'*reshape(h,k{i},n(i)*p(i)),n(i),p(i));
end

%% Calculate the auxiliary variables
for i=1:nperiods
  xa{i} = functions(i).fa(s{i},X{i},z{i},za{i},params);
  cxa{i} = funfitxy(fspace{i},Phi{i},xa{i});
end
if any(pa>0)
  %% With auxiliary expectations function
  if isfield(interp,'xa')
    xa = interp.xa;
    for i=1:nperiods, cxa{i} = funfitxy(fspace{i},Phi{i},xa{i}); end
  end
  xnrm     = inf;
  it       = 0;
  vec      = @(x) cell2mat(cellfun(@(z) z(:),x,'UniformOutput',false));
  
  [~,~,exitflag] = runeqsolver(@(X) ResidualVFI(X),xa{1}(:),...
                                -inf(numel(xa{1}),1),inf(numel(xa{1}),1),...
                                'krylov',struct('showiters',true,'Diagnostics'    , 'off' ,...
                              'DerivativeCheck', 'off' ,...
                                                'Jacobian'       , 'off',...
                                                'atol',1E-10,'rtol',1E-10));

  if showiters
    fprintf(1,'Solve for auxiliary variables\n');
    fprintf(1,'Successive approximation\n');
    fprintf(1,'   Iter\tResidual\n');
  end
  
  while(xnrm > atol && it < maxit)
    it    = it+1;
    xaold = xa;
    
    for i=nperiods:-1:1
      xanext = Phinext{i}*cxa{inext(i)};
      ha     = functions(i).ha(ss{i},xx{i},xa{i}(ind{i},:),ee{i},snext{i},xnext{i},xanext,params);
      za{i}  = reshape(shocks{i}.w'*reshape(ha,k{i},n(i)*pa(i)),n(i),pa(i));
      xa{i}  = functions(i).fa(s{i},X{i},z{i},za{i},params);
      cxa{i} = funfitxy(fspace{i},Phi{i},xa{i});
    end % for i=nperiods:-1:1
 
    xnrm = norm(vec(xa)-vec(xaold));
  if showiters && (it==1 || mod(it,showiters)==0 || xnrm <= atol)
    fprintf(1,'%7i\t%8.2E\n',it,xnrm); 
  end
 end % while(xnrm > atol && it < maxit)
end

%% Steady-state values
for i=1:nperiods
  xass{i} = funeval(cxa{i},fspace{i},model.ss.sss{i});
end

%% Export
interp.cxa  = cxa;
interp.xa   = xa;
interp.xass = xass;

%% Nested function
function R = ResidualVFI(xaold_first)
% RESIDUALVFI

  xaold_first = reshape(xaold_first,n(1),ma(1));
  cxa{1} = funfitxy(fspace{1},Phi{1},xaold_first);
  
  for i=nperiods:-1:1
    xanext = Phinext{i}*cxa{inext(i)};
    ha     = functions(i).ha(ss{i},xx{i},xa{i}(ind{i},:),ee{i},snext{i},xnext{i},xanext,params);
    za{i}  = reshape(shocks{i}.w'*reshape(ha,k{i},n(i)*pa(i)),n(i),pa(i));
    xa{i}  = functions(i).fa(s{i},X{i},z{i},za{i},params);
    cxa{i} = funfitxy(fspace{i},Phi{i},xa{i});
  end % for i=nperiods:-1:1
  
  
  R = xa{1}-xaold_first;
  R = R(:);

end

end