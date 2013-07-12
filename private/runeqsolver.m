function [x,f,exitflag] = runeqsolver(func,x,LB,UB,eqsolver,eqsolveroptions,varargin)
% RUNEQSOLVER Runs equations solvers (including MCP solvers)
%
% RUNEQSOLVER provides an unified interface to various equations and MCP
% solvers.
%
% X = RUNEQSOLVER(FUNC,X,LB,UB,EQSOLVER,EQSOLVEROPTIONS) tries to solve, using X
% as a starting point, the mixed complementarity problem of the form:
% LB =X     => FUNC(X)>0,
% LB<=X<=UB => FUNC(X)=0,
%     X =UB => FUNC(X)<0.
% LB and UB are the lower and upper bounds on X (not used if the chosen solver
% is not an MCP solver). RUNEQSOLVER returns X the solution. FUNC is an
% anonymous function that evaluates the equations, and possible the Jacobian at
% X.
% EQSOLVER designates the solver to use. Possible choices are 'lmmcp', 'fsolve',
% 'ncpsolve', and 'path'.
% EQSOLVEROPTIONS is a structure containing the options for the solver. See each
% solver's documentation for the available options. Common options are:
%    Jacobian        : 'off' to use numerical differentiation and 'on' if FUNC
%                      returns the Jacobian as second output.
%    DerivativeCheck : 'on' to check user-provided Jacobian against the one
%                      calculated by two-sided finite difference. 'off' to pass
%                      Jacobian check.
%
% X = RUNEQSOLVER(FUNC,X,LB,UB,EQSOLVER,EQSOLVEROPTIONS,VARARGIN) provides
% additional arguments for FUNC, which, in this case, takes the following form:
% FUNC(X,VARARGIN).
%
% [X,F] = RUNEQSOLVER(FUNC,X,LB,UB,EQSOLVER,EQSOLVEROPTIONS,...) returns F the
% value of the equations at solution.
%
% [X,F,EXITFLAG] = RUNEQSOLVER(FUNC,X,LB,UB,EQSOLVER,EQSOLVEROPTIONS,...)
% returns EXITFLAG, which describes the exit conditions. Possible values depend
% on the active solver, but a general rule is that 0 means failure to converge 1
% means convergence to a solution.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if strcmpi(eqsolver,'path'), global eqtosolve; end %#ok<TLEV>

if strcmpi(eqsolveroptions.Jacobian,'off') && ...
      any(strcmpi(eqsolver,{'lmmcp','ncpsolve','path'}))
  eqtosolve = @(Y) PbWithNumJac(func,Y,eqsolver,varargin);
else
  eqtosolve = @(Y) func(Y,varargin{:});
end


%% Derivative Check
if strcmpi(eqsolveroptions.DerivativeCheck,'on')
  Jnum  = numjac(func,x,struct([]),varargin{:});
  [~,J] = func(x,varargin{:});
  dJ    = norm(full(abs(J-Jnum)./max(1.0,abs(J))),Inf);
  fprintf(1,'Derivative Check: max relative diff = %g\n\n',dJ) %#ok<PRTCAL>
  eqsolveroptions.DerivativeCheck = 'off';
end


%% Solve equations
try
  switch eqsolver
    case 'lmmcp'
      [x,f,exitflag] = lmmcp(eqtosolve, x, LB, UB, eqsolveroptions);
    
    case 'sa'
      [x,f,exitflag] = SA(eqtosolve, x, eqsolveroptions);

    case 'krylov'
      [x,~,exitflag] = nsoli(eqtosolve, x, eqsolveroptions);
      exitflag = ~exitflag;
      f = NaN(size(x));

    case 'ncpsolve'
      exitflag = 1;  % ncpsolve does not output any exitflag on a failure
      [x,f]    = ncpsolve(@ncpsolvetransform, LB, UB, x, eqtosolve);
      f        = -f;
      if strcmp(lastwarn,'Failure to converge in ncpsolve')
        exitflag = 0;
        lastwarn('Warning reinitialization','RECS:WarningInit');
      end

    case 'path'
      [x,f,exitflag] = recspathmcp(x, LB, UB, 'pathtransform');
      clear global eqtosolve
    
    case 'mixed'
      reesolveroptions.maxit = 10;
      reesolveroptions.atol  = 1E-2;
      reesolveroptions.rtol  = 1E-3;
      [x,f,exitflag] = SA(eqtosolve, x, eqsolveroptions);
      
      reesolveroptions.maxit = 40;
      reesolveroptions.atol  = sqrt(eps);
      reesolveroptions.rtol  = sqrt(eps);
      [x,~,exitflag] = nsoli(eqtosolve, x, eqsolveroptions);
      exitflag = ~exitflag;
      f = NaN(size(x));

    case 'fsolve'
      options = optimset(optimset('Display','off'),eqsolveroptions);
      [x,f,exitflag] = fsolve(eqtosolve, x, options);

  end
catch
  exitflag = 0;
  f        = NaN(size(x));
end

function [F,J] = PbWithNumJac(func,Y,eqsolver,otherarg)
% PBWITHNUMJAC Calculates the numerical Jacobian of the problem func with respect to Y

if nargout==2
  J = numjac(func,Y,struct([]),otherarg{:});
  if strcmpi(eqsolver,'path'), J = sparse(J); end
end
F = func(Y,otherarg{:});
