function [x,f,exitflag] = runeqsolver(func,x,LB,UB,eqsolver,eqsolveroptions,varargin)
% RUNEQSOLVER Runs equations solvers (including MCP)
%
% RUNEQSOLVER is called by recsSolveDeterministicPb, recsSolveEquilibrium,
% recsSolveREEFull, recsSS. It is not meant to be called directly by the user.
  
% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if strcmp(eqsolver,'path'), global eqtosolve; end %#ok<TLEV>

if strcmp(eqsolveroptions.Jacobian,'off') && any(strcmp(eqsolver,{'lmmcp','ncpsolve','path'}))
  eqtosolve = @(Y) PbWithNumJac(func,Y,eqsolver,varargin);        
else
  eqtosolve = @(Y) func(Y,varargin{:});
end


%% Derivative Check
if strcmp(eqsolveroptions.DerivativeCheck,'on')
  Jnum  = numjac(func,x,[],varargin{:});
  [~,J] = func(x,varargin{:});
  dJ    = norm(full(abs(J-Jnum)./max(1.0,abs(J))),Inf);
  fprintf(1,'Derivative Check: max relative diff = %g\n\n',dJ) %#ok<PRTCAL>
  eqsolveroptions.DerivativeCheck = 'off';
end


%% Solve equations
switch eqsolver
  case 'lmmcp'
    [x,f,exitflag] = lmmcp(eqtosolve,...
                           x,LB,UB,...
                           eqsolveroptions);
    
  case 'fsolve'
    options = optimset(optimset('Display','off'),eqsolveroptions);
    [x,f,exitflag] = fsolve(eqtosolve,...
                            x,...
                            options);
 
  case 'ncpsolve'
    exitflag = 1;  % ncpsolve does not output any exitflag on a failure
    [x,f]    = ncpsolve(@ncpsolvetransform,...
                        LB,UB,x,...
                        eqtosolve);
    f        = -f;
  
  case 'path'
    [x,f,exitflag] = recspathmcp(x,LB,UB,'pathtransform');
    clear global eqtosolve

end

function [F,J] = PbWithNumJac(func,Y,eqsolver,otherarg)
% PBWITHNUMJAC Calculates the numerical Jacobian of the problem func with respect to Y

if nargout==2
  J = numjac(func,Y,[],otherarg{:});
  if strcmp(eqsolver,'path')
    J = sparse(J);
  end
end
F = func(Y,otherarg{:});
