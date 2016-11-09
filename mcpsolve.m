function [x,fval,exitflag,output,jacobian] = mcpsolve(f,x,lb,ub,options,varargin)
% MCPSOLVE solves a mixed complementarity problem
%
% MCPSOLVE is a reimplementation of the Fackler and Miranda's solver
% ncpsolve. MCPSOLVE adopts the standard convention for MCP problem, follows the
% convention of MATLAB optimization functions, and allows to solve array of
% problems.
%
% X = MCPSOLVE(F,X0) tries to solve the system of nonlinear equations F(X)=0 and
% starts at the vector X0. F accepts a vector X and return a vector F of equation
% values F evaluated at X and, as second output if required, a matrix J, the
% Jacobian evaluated at X.
%
% X = MCPSOLVE(F,X0,LB,UB) solves the mixed complementarity problem of the form:
% LB =X     =>   F(X)>0,
% LB<=X<=UB =>   F(X)=0,
%     X =UB =>   F(X)<0.
%
% [X,FVAL] = MCPSOLVE(F,X0,...) returns the value of the equations F at X.
%
% [X,FVAL,EXITFLAG] = MCPSOLVE(F,X0,...) returns EXITFLAG that describes the exit
% conditions. Possible values are
%      1         : MCPSOLVE converged to a root
%      0         : Too many iterations
%
% [X,FVAL,EXITFLAG,OUTPUT] = MCPSOLVE(F,X0,...) returns the structure OUTPUT
% that contains the number of function evaluations (OUTPUT.funcCount) and the
% number of iterations (OUTPUT.iterations).
%
% [X,FVAL,EXITFLAG,OUTPUT,JACOBIAN] = MCPSOLVE(F,X0,...) returns JACOBIAN the
% Jacobian of F evaluated at X.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
if nargin < 3 || isempty(lb), lb = -inf(size(x)); end
if nargin < 4 || isempty(ub), ub =  inf(size(x)); end

idlb  = isfinite(lb);
idub  = isfinite(ub);
ismcp = any([idlb(:); idub(:)]);

[nr,nc] = size(x);
n       = nr*nc;

defaultopt = struct(      ...
    'ArrayProblem'          , false,...
    'Display'               , 'none',...
    'MaxFunctionEvaluations', 100*nr*nc,...
    'MaxIterations'         , 400     ,...
    'maxsteps'              , 10,...
    'TolFun'                , 1E-10,...
    'TolX'                  , 1E-6,...
    'type'                  , 'smooth');
if nargin < 5 || isempty(options)
  options = defaultopt;
else
  options = catstruct(defaultopt,options);
end

ArrayProblem           = options.ArrayProblem;
MaxFunctionEvaluations = options.MaxFunctionEvaluations;
MaxIterations          = options.MaxIterations;
maxsteps               = max(0,options.maxsteps);
showiters              = strcmpi(options.Display,'iter');
TolFun                 = options.TolFun;
TolX                   = options.TolX;
if ~ismcp, type = 'notmcp';
else       type = options.type;
end

beta       = 0.5;
fnorm      = inf;
funcCount  = 0;
iterations = 0;
dx         = inf(size(x));
t          = 1;
if showiters
  fprintf(1,'                                   Norm of\n');
  fprintf(1,'Iteration Minor Func-count   f(x)    step\n');
  fprintf(1,'%5i      %2i %6i      %8.2E  %5.4g               (Input point)\n',0,0,0,[],[]);
end

while fnorm > TolFun && funcCount < MaxFunctionEvaluations && ...
      iterations < MaxIterations && norm(t*dx(:)) > TolX
  iterations = iterations + 1;
  xold = x;
  [fval,fjac] = feval(f,x,varargin{:});
  switch type
    case 'smooth'
      fval0 = fval;
      mcp_smooth_f();
    case 'minmax'
      mcp_minmax_f()
    case 'notmcp'
  end
  fnorm = norm(fval(:));
  if fnorm < TolFun, break, end
  if ~strcmp(type,'notmcp')
    if iterations==1
      if ~ArrayProblem
        if issparse(fjac), sputils = speye(n);
        else               sputils =   eye(n);
        end
      elseif nc~=1
        sputils = permute(repmat(eye(nc),1,1,nr),[3 1 2]);
      end
    end
    switch type
      case 'smooth'
        mcp_smooth_g();
      case 'minmax'
        mcp_minmax_g()
    end
  end
  if ArrayProblem
    if nc==1
      dx = -(fval./fjac);
    else
      dx = -arrayinv(fval,fjac);
    end
  else
    dx = -(fjac\fval);
  end
  x = xold + dx;
  fval = feval(f,x,varargin{:});
  switch type
    case 'smooth'
      mcp_smooth_f();
    case 'minmax'
      mcp_minmax_f()
    case 'notmcp'
  end
  fnormnew = norm(fval(:));

  %% Backstepping
  iterbackstep = 0;
  t = 1;
  while fnormnew > fnorm && iterbackstep < maxsteps
    iterbackstep = iterbackstep + 1;
    t = t*beta;
    x = xold + t*dx;
    fval = feval(f,x,varargin{:});
    switch type
      case 'smooth'
        mcp_smooth_f();
      case 'minmax'
        mcp_minmax_f()
      case 'notmcp'
    end
    fnormnew = norm(fval(:));
  end

  fnorm     = fnormnew;
  funcCount = funcCount + 2 + iterbackstep;
  if showiters
    fprintf(1,'%5i      %2i %6i      %8.2E  %5.4g\n',...
            [iterations iterbackstep funcCount fnorm t]);
  end
end

%% Output treatment
if nargout==5
  [fval,jacobian] = feval(f,x,varargin{:});
elseif nargout>=2
  fval = feval(f,x,varargin{:});
end

exitflag = 1;
output.funcCount  = funcCount;
output.iterations = iterations;

if iterations==MaxIterations
  exitflag = 0;
  if showiters
    fprintf(1,'Too many iterations\n');
  end
end

function mcp_smooth_f()

  % Phiplus
  fval(idub) = fval(idub) + (x(idub)-ub(idub)) + sqrt(fval(idub).^2+(x(idub)-ub(idub)).^2);
  % Phiminus
  fval(idlb) = fval(idlb) + (x(idlb)-lb(idlb)) - sqrt(fval(idlb).^2+(x(idlb)-lb(idlb)).^2);

end % mcp_smooth_f

function mcp_smooth_g()

  % Derivatives of phiplus
  sqplus = sqrt(fval0.^2+(x-ub).^2);

  dplus_du = 1 + fval0./sqplus;

  dplus_dv = zeros(size(x));
  dplus_dv(idub) = 1 + (x(idub)-ub(idub))./sqplus(idub);

  % Derivatives of phiminus
  phiplus = fval0;
  phiplus(idub) = phiplus(idub) + (x(idub)-ub(idub)) + sqrt(phiplus(idub).^2+(x(idub)-ub(idub)).^2);

  sqminus = sqrt(phiplus.^2+(x-lb).^2);

  dminus_du = 1-phiplus./sqminus;

  dminus_dv = zeros(size(x));
  dminus_dv(idlb) = 1 - (x(idlb)-lb(idlb))./sqminus(idlb);

  % Final computations
  if ArrayProblem
    if nc==1
      fjac = (dminus_du.*dplus_du).*fjac + dminus_du.*dplus_dv + dminus_dv;
    else
      inddiag            = 1:(nc+1):nc^2;
      sputils(:,inddiag) = dminus_du.*dplus_du;
      fjac               = arraymult(sputils,fjac,nr,nc,nc,nc);
      fjac(:,inddiag)    = fjac(:,inddiag) + (dminus_du.*dplus_dv+dminus_dv);
    end
  else
    inddiag = 1:(n+1):n^2;
    sputils(inddiag) = dminus_du.*dplus_du;
    fjac             = sputils*fjac;
    fjac(inddiag)    = fjac(inddiag) + (dminus_du.*dplus_dv+dminus_dv)';
  end % if ArrayProblem
end % mcp_smooth_g

function mcp_minmax_f()

  fval = min(max(fval,x-ub),x-lb);

end % mcp_minmax_f

function mcp_minmax_g()

  isatbounds = fval<(x-ub) | fval>(x-lb);
  if ArrayProblem 
    if nc==1
      fjac(isatbounds,:) = 1;
    else 
      sputils =   eye(nc);
      for i=1:nr
        fjac(i,isatbounds(i,:),:) = sputils(isatbounds(i,:),:);
      end
    end
  else
    fjac(isatbounds,:) = sputils(isatbounds,:);
  end

end % mcp_minmax_g

end