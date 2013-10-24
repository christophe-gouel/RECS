% TESTSJACOBIANS Tests the jacobian calculation in the variaous methods implemented in RECS

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

warning('off')

funapproxlist = {'resapprox' 'expfunapprox' 'expapprox'};


recsdirectory   = fileparts(which('recsSimul'));
addpath(fullfile(recsdirectory,'demos'))

expfunOK = [1 1 1 0 0];

for iter=1:5
  clear model interp s x
  switch iter
    case 1
      %% GRO1
      disp('Stochastic growth model')
      model = recsmodel('gro1.yaml',struct('Mu',0,'Sigma',0.007^2,'order',5));
      [interp,s] = recsinterpinit(10,[0.85*model.sss(1) min(model.shocks.e)*4],...
                                  [1.15*model.sss(1) max(model.shocks.e)*4],'cheb');
      [interp,x] = recsFirstGuess(interp,model,s,model.sss,model.xss);

    case 2
      %% GRO3
      disp('Stochastic growth model with recursive preferences and stochastic volatility')
      model = recsmodel('gro3.yaml',struct('Mu',[0 0],'Sigma',eye(2),'order',5));
      [interp,s] = recsinterpinit(4,...
                                  [0.85*model.sss(1) -0.11 log(0.007)*1.15],...
                                  [1.20*model.sss(1)  0.11 log(0.007)*0.85],...
                                  'cheb');
      [interp,x] = recsFirstGuess(interp,model,s);

    case 3
      %% CS1
      disp('Consumption/saving model with borrowing constraint')
      model      = recsmodel('cs1.yaml',struct('Mu',100,'Sigma',100,'order',5));
      [interp,s] = recsinterpinit(20,model.sss/2,model.sss*2);
      x          = s;

    case 4
      %% STO1
      disp('Competitive storage model')
      model = recsmodel('sto1.yaml',struct('Mu',1,'Sigma',0.05^2,'order',7));
      [interp,s] = recsinterpinit(40,model.sss*0.7,model.sss*1.5);
      [interp,x] = recsFirstGuess(interp,model,s,model.sss,model.xss,struct('T',5));

    case 5
      %% STO5
      disp('Two-country storage-trade model')
      model = recsmodel('sto5.yaml',struct('Mu',[1 1],...
                                           'Sigma',diag([0.05 0.05])^2,...
                                           'order',5));
      [interp,s] = recsinterpinit(7,0.73*model.sss,2*model.sss);
      [interp,x] = recsFirstGuess(interp,model,s,model.sss,model.xss,struct('T',5));

  end % switch

  [LB,UB] = model.functions.b(s,model.params);

  options = struct('display',0,...
                   'eqsolveroptions',struct('DerivativeCheck','on'),...
                   'reesolveroptions',struct('maxit',0));
  disp('Check equilibrium equations')
  options.reemethod = 'iter';
  for funapprox=funapproxlist
    if strcmp(funapprox,'expfunapprox') && ~expfunOK(iter), continue; end
    fprintf(1,' Functional approximation - %s\n',funapprox{1});
    options.funapprox = funapprox{1};
    recsSolveREE(interp,model,s,x,options);
  end

  disp('Check full Newton approach')
  options.reemethod = '1-step';
  for funapprox=funapproxlist
    if strcmp(funapprox,'expapprox'), continue; end
    if strcmp(funapprox,'expfunapprox') && ~expfunOK(iter), continue; end
    fprintf(1,' Functional approximation - %s\n',funapprox{1});
    options.funapprox = funapprox{1};
    recsSolveREE(interp,model,s,x,options);
  end

  if all(isinf([LB(:); UB(:)]))
    disp('Check iterative Newton approach')
    options.reesolver = 'fsolve';
    options.reemethod = 'iter';
    options.eqsolveroptions.DerivativeCheck = 'off';
    options.reesolveroptions.DerivativeCheck = 'on';
    options.reesolveroptions.MaxIter         = 0;
    options.reesolveroptions.Display         = 'off';
    for funapprox=funapproxlist
      if strcmp(funapprox,'expapprox'), continue; end
      if strcmp(funapprox,'expfunapprox') && ~expfunOK(iter), continue; end
      fprintf(1,' Functional approximation - %s\n',funapprox{1});
      options.funapprox = funapprox{1};
      recsSolveREE(interp,model,s,x,options);
    end
  end
end % iter

rmpath(fullfile(recsdirectory,'demos'))
