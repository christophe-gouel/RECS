warning('off')

reemethodlist = {'iter' 'iter-newton' '1-step'};
funapproxlist = {'resapprox-complete' 'resapprox-simple' 'expfunapprox' 'expapprox'};


recsdirectory   = fileparts(which('recsSimul'));
addpath(fullfile(recsdirectory,'demos'))
% CS1
% model      = recsmodelinit('cs1.yaml',struct('Mu',100,'Sigma',100,'order',5));
% [interp,s] = recsinterpinit(20,model.sss/2,model.sss*2);
% x          = s;

% GRO1
model = recsmodelinit('gro1.yaml',struct('Mu',0,'Sigma',0.007^2,'order',5));
[interp,s] = recsinterpinit(10,[0.85*model.sss(1) min(model.e)*4],...
                            [1.15*model.sss(1) max(model.e)*4],'cheb');
[interp,x] = recsFirstGuess(interp,model,s,model.sss,model.xss,50);

options = struct('display',0,...
                 'eqsolveroptions',struct('DerivativeCheck','on'),...
                 'reesolveroptions',struct('maxit',1));
disp('Check equilibrium equations')
options.reemethod = 'iter';
for funapprox=funapproxlist
  fprintf(1,' Functional approximation - %s\n',funapprox{1});
  options.funapprox = funapprox{1};
  recsSolveREE(interp,model,s,x,options);
end

disp('Check full Newton approach')
options.reemethod = '1-step';
for funapprox=funapproxlist
  fprintf(1,' Functional approximation - %s\n',funapprox{1});
  options.funapprox = funapprox{1};
  recsSolveREE(interp,model,s,x,options);
end

disp('Check iterative Newton approach')
options.reesolver = 'fsolve';
options.reemethod = 'iter-newton';
options.eqsolveroptions.DerivativeCheck = 'off';
options.reesolveroptions.DerivativeCheck = 'on';
options.reesolveroptions.MaxIter         = 0;
options.reesolveroptions.Display         = 'off';
for funapprox=funapproxlist
  fprintf(1,' Functional approximation - %s\n',funapprox{1});
  options.funapprox = funapprox{1};
  recsSolveREE(interp,model,s,x,options);
end

rmpath(fullfile(recsdirectory,'demos'))
