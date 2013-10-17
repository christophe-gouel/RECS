function [interp,X,exitflag,output] = recsSolveREESP(model,interp,X,options)
% RECSSOLVEREESP 

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                        ...
    'eqsolver'          , 'lmmcp'                          ,...
    'eqsolveroptions'   , struct('atol'           , sqrt(eps),...
                                 'Diagnostics'    , 'off'    ,...
                                 'DerivativeCheck', 'off'    ,...
                                 'Jacobian'       , 'on'     ,...
                                 'maxit'          , 1000)  ,...
    'extrapolate'       , 1                                ,...
    'Display'           , 'iter'                           ,...     
    'loop_over_s'       , 0                                ,...
    'funapprox'         , 'resapprox');
if nargin <=4
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(options,struct('funapprox','resapprox'));
  options = catstruct(defaultopt,options);
end

Display = options.Display;
switch lower(Display)
  case 'iter'
    showiters = 1;
  case {'off','none'}
    showiters = 0;
end
atol  = options.eqsolveroptions.atol;
maxit = options.eqsolveroptions.maxit;
extrapolate = options.extrapolate;

% Extract fields of model
nperiods  = model.nperiods;
params    = model.params;
shocks    = model.shocks;
functions = model.functions;
ixforward = cell(nperiods,1);
for i=1:nperiods, ixforward{i} = model.infos(i).ixforward; end
  
% Extract fields of interp
s          = interp.s;
fspace     = interp.fspace;
Phi        = interp.Phi;

cX = cell(nperiods,1);
for i=1:nperiods
  cX{i} = funfitxy(fspace{i},Phi{i},X{i});
end

cnrm     = 1;
it       = 0;
vec      = @(x) cell2mat(cellfun(@(z) reshape(z,[],1),x,'UniformOutput',false));
inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);
exitEQ   = zeros(nperiods,1);

%% Solve for the rational expectations equilibrium
if showiters
  fprintf(1,'Successive approximation\n');
  fprintf(1,'   Iter\tResidual\n');
end

while(cnrm > atol && it < maxit)
  it    = it+1;
  cXold = cX;
  
  for i=nperiods:-1:1
    [LB,UB] = functions(i).b(s{i},params);
    [X{i},~,exitEQ(i)] = recsSolveEquilibrium(s{i},X{i},[],...
                                             functions(inext(i)).b,...
                                             functions(i).f,...
                                             functions(i).g,...
                                             functions(i).h,...
                                             params,...
                                             cX{inext(i)}(:,ixforward{i}),...
                                             shocks{i}.e,shocks{i}.w,...
                                             fspace{inext(i)},...
                                             ixforward{i},...
                                             options,LB,UB);

    cX{i} = funfitxy(fspace{i},Phi{i},X{i});
  end % for i=nperiods:-1:1
  cnrm = norm(vec(cX)-vec(cXold));
  if showiters, fprintf(1,'%7i\t%8.2E\n',it,cnrm); end
end % while(cnrm > atol && it < maxit)

%% Output treatment
% Check if state satisfies bounds
output = struct();
output.snextmin = cell(nperiods,1);
output.snextmax = cell(nperiods,1);
for i=1:nperiods
  si     = s{i};
  xi     = X{i};
  ei     = shocks{i}.e;
  n      = size(si,1);
  k      = size(ei,1);
  ind    = (1:n);
  ind    = ind(ones(1,k),:);
  ss     = si(ind,:);
  ee     = ei(repmat(1:k,1,n),:);
  xx     = xi(ind,:);
  snext  = functions(i).g(ss,xx,ee,params);
  output.snextmin{i} = min(snext);
  output.snextmax{i} = max(snext);
  if extrapolate==2 || extrapolate==-1
    vari   = 1:fspace{inext(i)}.d;
    varmin = vari(min(snext)<fspace{inext(i)}.a);
    varmax = vari(max(snext)>fspace{inext(i)}.b);
    if ~isempty(varmin)
      warning('RECS:Extrapolation','State variables (%s) in subperiod %i beyond smin',...
              int2str(varmin),inext(i))
    end
    if ~isempty(varmax)
      warning('RECS:Extrapolation','State variables (%s) in subperiod %i beyond smax',...
              int2str(varmax),inext(i))
    end
  end
end % for i=1:nperiods

% exitflag
if it==maxit
  exitflag = 0;
  if showiters
    warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
    fprintf(1,'Too many iterations\n'); 
  end
else
  if all(exitEQ)
    exitflag = 1;
    if showiters
      fprintf(1,'Solution found - Residual lower than absolute tolerance\n');
    end
  else
    exitflag = 0;
    warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
  end
end

% interp
interp.cX = cX;
interp.X  =  X;

