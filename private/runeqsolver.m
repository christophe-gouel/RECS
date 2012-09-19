function [x,f,exitflag] = runeqsolver(func,x,LB,UB,eqsolver,eqsolveroptions,varargin)
% RUNEQSOLVER Runs equations solvers (including MCP)
%
% RUNEQSOLVER is called by recsSolveDeterministicPb, recsSolveEquilibrium,
% recsSolveREEFull, recsSS. It is not meant to be called directly by the user.
  
% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Derivative Check
if strcmp(eqsolveroptions.DerivativeCheck,'on')
  Jnum  = numjac(func,x,[],varargin{:});
  [~,J] = func(x,varargin{:});
  dJ    = norm(full(abs(J-Jnum)./max(1.0,abs(J))),Inf);
  fprintf(1,'Derivative Check: max relative diff = %g\n\n',dJ)
  eqsolveroptions.DerivativeCheck = 'off';
end

%% Solve equations
switch eqsolver
  case 'fsolve'
    options = optimset('Display','off',...
                       'Jacobian','on');
    options = optimset(options,eqsolveroptions);
    [x,f,exitflag] = fsolve(@(Y) func(Y,varargin{:}),...
                            x,...
                            options);
 
  case 'lmmcp'
    [x,f,exitflag] = lmmcp(@(Y) func(Y,varargin{:}),...
                           x,...
                           LB,...
                           UB,...
                           eqsolveroptions);
 
  case 'ncpsolve'
    exitflag = 1;  % ncpsolve does not output any exitflag on a failure
    [x,f]    = ncpsolve(@ncpsolvetransform,...
                        LB,...
                        UB,...
                        x,...
                        func,...
                        varargin{:});
    f        = -f;
 
  case 'path'
    global par %#ok<TLEV>
    par        = [{func} varargin];
    [x,f,exitflag] = pathmcp(x,LB,UB,'pathtransform');
    clear global par

end
