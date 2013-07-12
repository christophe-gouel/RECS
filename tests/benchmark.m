warning('off','backtrace')
addpath('../demos');
% cs1
mc1 = recsmodel('cs1.yaml',struct('Mu',100,'Sigma',10^2,'order',5));
[ic1,sc1] = recsinterpinit(20,mc1.sss/2,mc1.sss*2);
Tc1 = 15;
% sto1
ms1 = recsmodel('sto1.yaml',struct('Mu',1,'Sigma',0.05^2,'order',7));
[is1,ss1] = recsinterpinit(30,ms1.sss*0.7,ms1.sss*1.5);
Ts1 = 5;
% sto2
ms2 = recsmodel('sto2.yaml',struct('Mu',1,'Sigma',0.05^2,'order',7));
[is2,ss2] = recsinterpinit(50,0.7,2);
Ts2 = 5;
% gro1
mg1 = recsmodel('gro1.yaml',struct('Mu',0,'Sigma',0.007^2,'order',5));
[ig1,sg1] = recsinterpinit(10,[0.85*mg1.sss(1) min(mg1.e)*4],[1.15*mg1.sss(1) max(mg1.e)*4]);
Tg1 = 50;
% gro2
mg2 = recsmodel('gro2.yaml',struct('Mu',0,'Sigma',0.04^2,'order',7));
[ig2,sg2] = recsinterpinit(14,[0.47*mg2.sss(1)  min(mg2.e)*3.5],[1.72*mg2.sss(1)  max(mg2.e)*3.5]);
Tg2 = 50;

mcpsolverlist    = {'lmmcp','ncpsolve','path'};
reemethodlist    = {'iter-newton','1-step'};
modellist        = {'c1','g1','g2','s1','s2'};
funapproxlist    = {'resapprox','expapprox'};

time = zeros(length(modellist)*length(mcpsolverlist),4);

solvervar = cell(length(modellist)*length(mcpsolverlist),1);
modelvar  = cell(length(modellist)*length(mcpsolverlist),1);

for j=1:2
  i = 0;
  for modeliter=modellist
    interp = eval(['i' modeliter{1}]);
    model  = eval(['m' modeliter{1}]);
    T      = eval(['T' modeliter{1}]);
    s      = eval(['s' modeliter{1}]);
    for eqsolver=mcpsolverlist
      i = i+1;
      solvervar{i} = eqsolver{1};
      modelvar{i}  = modeliter{1};
      options.eqsolver  = eqsolver{1};
      options.funapprox = funapproxlist{1};
      tic; [interp,~,~,exitflag] = recsFirstGuess(interp,model,s,model.sss,model.xss,Ts1,options); tmp = toc;
      if exitflag==1, time(i,1) = tmp; else time(i,1) = NaN; end
      options.reemethod = reemethodlist{1};
      tic; [~,~,~,~,exitflag] = recsSolveREE(interp,model,s,[],options); tmp = toc;
      if exitflag==1, time(i,2) = tmp; else time(i,2) = NaN; end
      options.reemethod = reemethodlist{2};
      tic; [~,~,~,~,exitflag] =  recsSolveREE(interp,model,s,[],options); tmp = toc;
      if exitflag==1, time(i,3) = tmp; else time(i,3) = NaN; end
      options.reemethod = reemethodlist{1};
      options.funapprox = funapproxlist{2};
      tic; [~,~,~,~,exitflag] = recsSolveREE(interp,model,s,[],options); tmp = toc;
      if exitflag==1, time(i,4) = tmp; else time(i,4) = NaN; end
    end
  end
end

Res = dataset(modelvar,solvervar,time(:,1),time(:,2),time(:,3),time(:,4));
Res.Properties.VarNames = {'Model',...
                           'Solver',...
                           'First_guess',...
                           'Response_approx_Iter', ...
                           'Response_approx_1step',...
                           'Expectations_approx'};
