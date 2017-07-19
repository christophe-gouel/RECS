function [se,lEf,lEf_res,lEE,lEE_res] = recsAccuracySP(model,interp,se,options)
% RECSACCURACYSP Evaluates approximation accuracy

%% Initialization
% Options
defaultopt = struct(...
    'display'         , 1                                  ,...
    'eqsolver'        , 'lmmcp'                            ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off'    ,...
                               'DerivativeCheck', 'off'    ,...
                               'Jacobian'       , 'on')    ,...
    'loop_over_s'     , 0                                  ,...
    'simulmethod'     , 'interpolation'                    ,...
    'UseParallel'     , 'always');
if nargin<4
  options = defaultopt;
else
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end
display     = options.display;

functions = model.functions;
shocks    = model.shocks;
params    = model.params;
nperiods  = model.nperiods;
ixforward = cell(nperiods,1);
for i=1:nperiods, ixforward{i} = model.infos(i).ixforward; end

cX     = interp.cX;
fspace = interp.fspace;

inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);

n      = size(se{1},1);

xe         = cell(nperiods,1);
for i=1:nperiods
  [LB,UB] = functions(i).b(se{i},params);
  xe{i}   = min(max(funeval(cX{i},fspace{i},se{i}),LB),UB);

  if strcmp(options.simulmethod,'solve')
    xe{i} = recsSolveEquilibrium(se{i},...
                                 xe{i},...
                                 zeros(size(xe{i},1),0),...
                                 functions(inext(i)).b,...
                                 functions(i).f,...
                                 functions(i).g,...
                                 functions(i).h,...
                                 params,...
                                 cX{inext(i)}(:,ixforward{i}),...
                                 shocks{i}.e,shocks{i}.w,...
                                 fspace{inext(i)},...
                                 ixforward{i},options,LB,UB);
  end
end

%% Calculation of ze
ze         = cell(nperiods,1);
for i=1:nperiods
  e         = shocks{i}.e;
  k         = length(e);
  p         = model.dim{i,3};
  ind       = (1:n);
  ind       = ind(ones(1,k),:);
  ss        = se{i}(ind,:);
  xx        = xe{i}(ind,:);
  ee        = e(repmat(1:k,1,n),:);
  sen       = functions(i).g(ss,xx,ee,params);
  if i~=nperiods, j = i+1;
  else,           j = 1;
  end
  [LBn,UBn] = functions(j).b(sen,params);
  xen       = min(max(funeval(cX{j},fspace{j},sen),LBn),UBn);
  hv        = functions(i).h(ss,xx,ee,sen,xen,params);
  ze{i}     = reshape(shocks{i}.w'*reshape(hv,k,n*p),n,p);
end

%% Equilibrium equation error
Ef = cell(nperiods,1); lEf = Ef; lEf_res = Ef;
for i=1:nperiods
  fe        = functions(i).f(se{i},xe{i},ze{i},params);
  [LB,UB]   = functions(i).b(se{i},params);
  Ef{i}     = abs(min(max(-fe,LB-xe{i}),UB-xe{i}));
  lEf{i}    = log10(Ef{i});
  lEf_res{i} = [log10(max(Ef{i}));
                log10(sum(Ef{i})/size(Ef{i},1))];
end

%% Euler equation error
EE = cell(nperiods,1); lEE = EE; lEE_res = EE;
for i=1:nperiods
  EE{i}      = functions(i).ee(se{i},xe{i},ze{i},params);
  lEE{i}     = log10(abs(EE{i}));
  lEE_res{i} = [log10(max(abs(EE{i})));
                log10(sum(abs(EE{i}))/size(EE{i},1))];
end

% Works only if there is the output of functions.ee is unidimensional
if display==1 && all(~isnan(cell2mat(lEE_res)))
  lEE_res = cell2mat(lEE_res);
  lEE_res = reshape(lEE_res,[],nperiods)';
  lEE_res = array2table(lEE_res,'VariableNames',{'Max' 'Mean'});
  lEE_res.Season = (1:nperiods)';
  disp(' Euler equation error (in log10)');
  disp(lEE_res);
end
