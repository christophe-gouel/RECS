function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = lmmcp(FUN,x,lb,ub,options,varargin)
% LMMCP solves a mixed complementarity problem.
%
% LMMCP uses a semismooth least squares formulation. The method applies a
% Levenberg-Marquardt/Gauss-Newton algorithm to a least-squares formulation.
%
% X = LMMCP(FUN,X0) tries to solve the system of nonlinear equations F(X)=0 and
% starts at the vector X0. FUN accepts a vector X and return a vector F of equation
% values F evaluated at X and, as second output if required, a matrix J, the
% Jacobian evaluated at X.
%
% X = LMMCP(FUN,X0,LB,UB) solves the mixed complementarity problem of the form:
% LB =X     =>   F(X)>0,
% LB<=X<=UB =>   F(X)=0,
%     X =UB =>   F(X)<0.
%
% X = LMMCP(FUN,X0,LB,UB,OPTIONS) solves the MCP problem using the options
% defined in the structure OPTIONS. Main fields are
%      Display    : control the display of iterations, 'none' (default),
%                   'iter-detailed' or 'final-detailed'
%  Switch from phase I to phase II
%      preprocess : activate preprocessor for phase I (default = 1)
%      presteps   : number of iterations in phase I (default = 20)
%  Termination parameters
%      MaxIter    : Maximum number of iterations (default = 500)
%      tmin       : safeguard stepsize (default = 1E-12)
%      TolFun     : Termination tolerance on the function value, a positive 
%                   scalar (default = sqrt(eps))
%  Stepsize parameters
%      m          : number of previous function values to use in the nonmonotone
%                   line search rule (default = 10)
%      kwatch     : maximum number of steps (default = 20 and should not be
%                   smaller than m)
%      watchdog   : activate the watchdog strategy (default = 1)
%  Ther are other minor parameters. Please see the code for their default values
%  and interpretation.
%
% [X,FVAL] = LMMCP(FUN,X0,...) returns the value of the equations FUN at X.
%
% [X,FVAL,EXITFLAG] = LMMCP(FUN,X0,...) returns EXITFLAG that describes the exit
% conditions. Possible values are
%      1         : LMMCP converged to a root
%      0         : Too many iterations
%     -1         :
%
% [X,FVAL,EXITFLAG,OUTPUT] = LMMCP(FUN,X0,...) returns the structure OUTPUT that
% contains the number of iterations (OUTPUT.iterations), the value of the merit
% function (OUTPUT.Psix), and the norm of the derivative of the merit function
% (OUTPUT.normDPsix).
%
% [X,FVAL,EXITFLAG,OUTPUT,JACOB] = LMMCP(FUN,X0,...) returns JACOB the Jacobian
% of FUN evaluated at X.
%
% More details of the main program may be found in the following paper:
%
% Christian Kanzow and Stefania Petra: On a semismooth least squares formulation of
% complementarity problems with gap reduction. Optimization Methods and Software
% 19, 2004, pp. 507-525.
%
% In addition, the current implementation uses a preprocessor which is the
% projected Levenberg-Marquardt step from the following preprint:
%
% Christian Kanzow and Stefania Petra: Projected filter trust region methods for a
% semismooth least squares formulation of mixed complementarity
% problems. Optimization Methods and Software
% 22, 2007, pp. 713-735.
%
% A user's guide is also available:
%
% Christian Kanzow and Stefania Petra (2005).
% LMMCP --- A Levenberg-Marquardt-type MATLAB Solver for Mixed Complementarity Problems.
% University of Wuerzburg.
% http://www.mathematik.uni-wuerzburg.de/~kanzow/software/UserGuide.pdf
%
% This is a modification by Christophe Gouel of the original files, which can be
% downloaded from:
% http://www.mathematik.uni-wuerzburg.de/~kanzow/software/LMMCP.zip
%
% copyright: Christian Kanzow and Stefania Petra
%            Institute of Applied Mathematics and Statistics
%            University of Wuerzburg
%            Am Hubland
%            97074 Wuerzburg
%            GERMANY
%
%            e-mail: kanzow@mathematik.uni-wuerzburg.de
%                    petra@mathematik.uni-wuerzburg.de

%% Initialization
defaultopt = struct(...
    'beta',       0.55,...
    'Big',        1e10,...
    'delta',      5,...
    'deltamin',   1,...
    'Display',    'none',...
    'epsilon1',   1e-6,...
    'eta',        0.95,...
    'kwatch',     20,...
    'lambda1',    0.1,...
    'm',          10,...
    'MaxIter',    500,...
    'null',       1e-10,...
    'preprocess', 1,...
    'presteps',   20,...
    'sigma',      1e-4,...
    'sigma1',     0.5,...
    'sigma2',     2,...
    'tmin',       1e-12,...
    'TolFun',     sqrt(eps),...
    'watchdog',   1);

if nargin < 4
  ub = inf(size(x));
  if nargin < 3
    lb = -inf(size(x));
  end
end

if nargin < 5 || isempty(options) || ~isstruct(options)
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

warning('off','MATLAB:rankDeficientMatrix')

switch options.Display
  case {'off','none'}
    verbosity = 0;
  case {'iter','iter-detailed'}
    verbosity = 2;
  case {'final','final-detailed'}
    verbosity = 1;
  otherwise
    verbosity = 0;
end

% parameter settings
eps1 = options.epsilon1;
eps2 = 0.5*options.TolFun^2;
null = options.null;
Big  = options.Big;

% maximal number of iterations
kmax     = options.MaxIter;

% choice of lambda
lambda1  = options.lambda1;
lambda2  = 1-lambda1;

% steplength parameters
beta     = options.beta;
sigma    = options.sigma;
tmin     = options.tmin;

% parameters watchdog and nonmonotone line search; redefined later
m        = options.m;
kwatch   = options.kwatch;
watchdog = options.watchdog; % 1=watchdog strategy active, otherwise not

% parameters for preprocessor
preprocess = options.preprocess; % 1=preprocessor used, otherwise not
presteps   = options.presteps; % maximum number of preprocessing steps

% trust-region parameters for preprocessor
delta    = options.delta;
deltamin = options.deltamin;
sigma1   = options.sigma1;
sigma2   = options.sigma2;
eta      = options.eta;

% initializations
k        = 0;
k_main   = 0;

% compute a feasible starting point by projection
x = max(lb,min(x,ub));

n = length(x);
OUTPUT.Dim = n;

% definition of index sets I_l, I_u and I_lu
Indexset       = zeros(n,1);
I_l            = lb>-Big & ub>Big;
I_u            = lb<-Big & ub<Big;
I_lu           = lb>-Big & ub<Big;
Indexset(I_l)  = 1;
Indexset(I_u)  = 2;
Indexset(I_lu) = 3;

% function evaluations
[Fx,DFx] = feval(FUN,x,varargin{:});

% choice of NCP-function and corresponding evaluations
Phix      = Phi(x,Fx,lb,ub,lambda1,lambda2,n,Indexset);
normPhix  = norm(Phix);
Psix      = 0.5*(Phix'*Phix);
DPhix     = DPhi(x,Fx,DFx,lb,ub,lambda1,lambda2,n,Indexset);
DPsix     = DPhix'*Phix;
normDPsix = norm(DPsix);

% save initial values
x0         = x;
Phix0      = Phix;
Psix0      = Psix;
DPhix0     = DPhix;
DPsix0     = DPsix;
normDPsix0 = normDPsix;

% watchdog strategy
aux    = zeros(m,1);
aux(1) = Psix;
MaxPsi = Psix;

if watchdog==1
  kbest        = k;
  xbest        = x;
  Phibest      = Phix;
  Psibest      = Psix;
  DPhibest     = DPhix;
  DPsibest     = DPsix;
  normDPsibest = normDPsix;
end

% initial output
if verbosity > 1
  fprintf('   k               Psi(x)                || DPsi(x) ||    stepsize\n');
  disp('====================================================================')
  disp('********************* Output at starting point *********************')
  fprintf('%4.0f %24.5e %24.5e\n',k,Psix,normDPsix);
end

%% Preprocessor using local method

if preprocess==1

  if verbosity > 1
    disp('************************** Preprocessor ****************************')
  end
  
  normpLM=1;
  while (k < presteps) && (Psix > eps2) && (normpLM>null)
    k = k+1;
    
    % choice of Levenberg-Marquardt parameter, note that we do not use
    % the condition estimator for large-scale problems, although this
    % may cause numerical problems in some examples

    i  = false;
    mu = 0;
    if n<100
      i = true;
      mu = 1e-16;
      if condest(DPhix'*DPhix)>1e25
        mu = 1e-6/(k+1);
      end
    end
    if i
      pLM =  [DPhix; sqrt(mu)*speye(n)]\[-Phix; zeros(n,1)];
    else
      pLM = -DPhix\Phix;
    end
    normpLM = norm(pLM);
    
    % compute the projected Levenberg-Marquard step onto box Xk
    lbnew = max(min(lb-x,0),-delta);
    ubnew = min(max(ub-x,0),delta);
    d     = max(lbnew,min(pLM,ubnew));
    xnew  = x+d;

    % function evaluations etc.
    [Fxnew,DFxnew] = feval(FUN,xnew,varargin{:});
    Phixnew        = Phi(xnew,Fxnew,lb,ub,lambda1,lambda2,n,Indexset);
    Psixnew        = 0.5*(Phixnew'*Phixnew);
    normPhixnew    = norm(Phixnew);
    
    % update of delta
    if normPhixnew<=eta*normPhix
      delta = max(deltamin,sigma2*delta);
    elseif normPhixnew>5*eta*normPhix
      delta = max(deltamin,sigma1*delta);
    end

    % update
    x         = xnew;
    Fx        = Fxnew;
    DFx       = DFxnew;
    Phix      = Phixnew;
    Psix      = Psixnew;
    normPhix  = normPhixnew;
    DPhix     = DPhi(x,Fx,DFx,lb,ub,lambda1,lambda2,n,Indexset);
    DPsix     = DPhix'*Phix;
    normDPsix = norm(DPsix,inf);
    
    % output at each iteration
    t=1;
    if verbosity > 1
      fprintf('%4.0f %24.5e %24.5e %11.7g\n',k,Psix,normDPsix,t);
    end
  end
end

% terminate program or redefine current iterate as original initial point
if preprocess==1 && Psix<eps2
  if verbosity > 0
    fprintf('Psix = %1.4e\nnormDPsix = %1.4e\n',Psix,normDPsix);
    disp('Approximate solution found.')
  end
  EXITFLAG          = 1;
  FVAL              = Fx;
  OUTPUT.iterations = k;
  OUTPUT.Psix       = Psix;
  OUTPUT.normDPsix  = normDPsix;
  JACOB             = DFx;
  return
elseif preprocess==1 && Psix>=eps2
  x=x0;
  Phix=Phix0;
  Psix=Psix0;
  DPhix=DPhix0;
  DPsix=DPsix0;
  if verbosity > 1
    disp('******************** Restart with initial point ********************')
    fprintf('%4.0f %24.5e %24.5e\n',k_main,Psix0,normDPsix0);
  end
end

%%   Main algorithm

if verbosity > 1
  disp('************************** Main program ****************************')
end

while (k < kmax) && (Psix > eps2)

  % choice of Levenberg-Marquardt parameter, note that we do not use
  % the condition estimator for large-scale problems, although this
  % may cause numerical problems in some examples

  i = false;
  if n<100
    i  = true;
    mu = 1e-16;
    if condest(DPhix'*DPhix)>1e25
      mu = 1e-1/(k+1);
    end
  end
  
  % compute a Levenberg-Marquard direction

  if i
    d = [DPhix; sqrt(mu)*speye(n)]\[-Phix; zeros(n,1)];
  else
    d = -DPhix\Phix;
  end

  % computation of steplength t using the nonmonotone Armijo-rule
  % starting with the 6-th iteration

  % computation of steplength t using the monotone Armijo-rule if
  % d is a 'good' descent direction or k<=5

  t       = 1;
  xnew    = x+d;
  Fxnew   = feval(FUN,xnew,varargin{:});
  Phixnew = Phi(xnew,Fxnew,lb,ub,lambda1,lambda2,n,Indexset);
  Psixnew = 0.5*(Phixnew'*Phixnew);
  const   = sigma*DPsix'*d;
  
  while (Psixnew > MaxPsi + const*t)  && (t > tmin)
    t       = t*beta;
    xnew    = x+t*d;
    Fxnew   = feval(FUN,xnew,varargin{:});
    Phixnew = Phi(xnew,Fxnew,lb,ub,lambda1,lambda2,n,Indexset);
    Psixnew = 0.5*(Phixnew'*Phixnew);
  end

  % updatings
  x         = xnew;
  Fx        = Fxnew;
  Phix      = Phixnew;
  Psix      = Psixnew;
  [~,DFx]   = feval(FUN,x,varargin{:});
  DPhix     = DPhi(x,Fx,DFx,lb,ub,lambda1,lambda2,n,Indexset);
  DPsix     = DPhix'*Phix;
  normDPsix = norm(DPsix);
  k         = k+1;
  k_main    = k_main+1;
  
  if k_main<=5
    aux(mod(k_main,m)+1) = Psix;
    MaxPsi               = Psix;
  else
    aux(mod(k_main,m)+1) = Psix;
    MaxPsi               = max(aux);
  end
  
  % updatings for the watchdog strategy
  if watchdog ==1
    if Psix<Psibest
      kbest        = k;
      xbest        = x;
      Phibest      = Phix;
      Psibest      = Psix;
      DPhibest     = DPhix;
      DPsibest     = DPsix;
      normDPsibest = normDPsix;
    elseif k-kbest>kwatch
      x=xbest;
      Phix=Phibest;
      Psix=Psibest;
      DPhix=DPhibest;
      DPsix=DPsibest;
      normDPsix=normDPsibest;
      MaxPsi=Psix;
    end
  end

  if verbosity > 1
    % output at each iteration
    fprintf('%4.0f %24.5e %24.5e %11.7g\n',k,Psix,normDPsix,t);
  end
end

%% Final output
if Psix<=eps2
  EXITFLAG = 1;
  if verbosity > 0, disp('Approximate solution found.'); end
elseif k>=kmax
  EXITFLAG = 0;
  if verbosity > 0, disp('Maximum iteration number reached.'); end
elseif normDPsix<=eps1
  EXITFLAG          = -1; % Provisoire
  if verbosity > 0, disp('Approximate stationary point found.'); end
else
  EXITFLAG          = -1; % Provisoire
  if verbosity > 0, disp('No solution found.'); end
end

FVAL              = Fx;
OUTPUT.iterations = k;
OUTPUT.Psix       = Psix;
OUTPUT.normDPsix  = normDPsix;
JACOB             = DFx;

%% Subfunctions

function y = Phi(x,Fx,lb,ub,lambda1,lambda2,n,Indexset)
%% PHI

y           = zeros(2*n,1);
phi_u       = sqrt((ub-x).^2+Fx.^2)-ub+x+Fx;
LZ          = false(n,1); % logical zero

I0          = Indexset==0;
y(I0)       = -lambda1*Fx(I0);
y([LZ; I0]) = -lambda2*Fx(I0);

I1          = Indexset==1;
y(I1)       = lambda1*(-x(I1)+lb(I1)-Fx(I1)+sqrt((x(I1)-lb(I1)).^2+Fx(I1).^2));
y([LZ; I1]) = lambda2*max(0,x(I1)-lb(I1)).*max(0,Fx(I1));

I2          = Indexset==2;
y(I2)       = -lambda1*phi_u(I2);
y([LZ; I2]) = lambda2*max(0,ub(I2)-x(I2)).*max(0,-Fx(I2));

I3          = Indexset==3;
y(I3)       = lambda1*(sqrt((x(I3)-lb(I3)).^2+phi_u(I3).^2)-x(I3)+lb(I3)-phi_u(I3));
y([LZ; I3]) = lambda2*(max(0,x(I3)-lb(I3)).*max(0,Fx(I3))+max(0,ub(I3)-x(I3)).*max(0,-Fx(I3)));


function H = DPhi(x,Fx,DFx,lb,ub,lambda1,lambda2,n,Indexset)
%% DPHI evaluates an element of the C-subdifferential of operator Phi

null       = 1e-8;
beta_l     = zeros(n,1);
beta_u     = zeros(n,1);
alpha_l    = zeros(n,1);
alpha_u    = zeros(n,1);


z          = zeros(n,1);
H2         = sparse(n,n);

I          = abs(x-lb)<=null & abs(Fx)<=null;
beta_l(I)  = 1;
z(I)       = 1;

I          = abs(ub-x)<=null & abs(Fx)<=null;
beta_u(I)  = 1;
z(I)       = 1;

I          = x-lb>=-null & Fx>=-null;
alpha_l(I) = 1;

I          = ub-x>=-null & Fx<=null;
alpha_u(I) = 1;

Da         = zeros(n,1);
Db         = zeros(n,1);

I          = 1:n;

I0         = Indexset==0;
Da(I0)     = 0;
Db(I0)     = -1;
H2(I0,:)   = -DFx(I0,:);

I1         = Indexset==1;
denom1     = zeros(n,1);
denom2     = zeros(n,1);
if any(I1)
  denom1(I1) = max(null,sqrt((x(I1)-lb(I1)).^2+Fx(I1).^2));
  denom2(I1) = max(null,sqrt(z(I1).^2+(DFx(I1,:)*z).^2));
end

I1b        = Indexset==1 & beta_l==0;
Da(I1b)    = (x(I1b)-lb(I1b))./denom1(I1b)-1;
Db(I1b)    = Fx(I1b)./denom1(I1b)-1;
I1b        = Indexset==1 & beta_l~=0;
if any(I1b)
  Da(I1b)  = z(I1b)./denom2(I1b)-1;
  Db(I1b)  = (DFx(I1b,:)*z)./denom2(I1b)-1;
end

I1a         = I(Indexset==1 & alpha_l==1);
if any(I1a)
  H2(I1a,:) = bsxfun(@times,x(I1a)-lb(I1a),DFx(I1a,:))+...
              sparse(1:length(I1a),I1a,Fx(I1a),length(I1a),n,length(I1a));
end

I2         = Indexset==2;
denom1     = zeros(n,1);
denom2     = zeros(n,1);
if any(I2)
  denom1(I2) = max(null,sqrt((ub(I2)-x(I2)).^2+Fx(I2).^2));
  denom2(I2) = max(null,sqrt(z(I2).^2+(DFx(I2,:)*z).^2));
end

I2b        = Indexset==2 & beta_u==0;
Da(I2b)    = (ub(I2b)-x(I2b))./denom1(I2b)-1;
Db(I2b)    = -Fx(I2b)./denom1(I2b)-1;
I2b        = Indexset==2 & beta_u~=0;
if any(I2b)
  Da(I2b)  = -z(I2b)./denom2(I2b)-1;
  Db(I2b)  = -(DFx(I2b,:)*z)./denom2(I2b)-1;
end

I2a         = I(Indexset==2 & alpha_u==1);
if any(I2a)
  H2(I2a,:) = bsxfun(@times,x(I2a)-ub(I2a),DFx(I2a,:))+...
              sparse(1:length(I2a),I2a,Fx(I2a),length(I2a),n,length(I2a));
end

I3         = Indexset==3;
phi        = zeros(n,1);
ai         = zeros(n,1);
bi         = zeros(n,1);
ci         = zeros(n,1);
di         = zeros(n,1);
denom1     = zeros(n,1);
denom2     = zeros(n,1);
denom3     = zeros(n,1);
denom4     = zeros(n,1);
if any(I3)
  phi(I3)    = -ub(I3)+x(I3)+Fx(I3)+sqrt((ub(I3)-x(I3)).^2+Fx(I3).^2);
  denom1(I3) = max(null,sqrt((x(I3)-lb(I3)).^2+phi(I3).^2));
  denom2(I3) = max(null,sqrt(z(I3).^2+(DFx(I3,:)*z).^2));
  denom3(I3) = max(null,sqrt((ub(I3)-x(I3)).^2+Fx(I3).^2));
  denom4(I3) = max(null,sqrt(z(I3).^2));
end

I3bu       = Indexset==3 & beta_u==0;
ci(I3bu)   = (x(I3bu)-ub(I3bu))./denom3(I3bu)+1;
di(I3bu)   = Fx(I3bu)./denom3(I3bu)+1;
I3bu       = Indexset==3 & beta_u~=0;
if any(I3bu)
  ci(I3bu)   = 1+z(I3bu)./denom2(I3bu);
  di(I3bu)   = 1+(DFx(I3bu,:)*z)./denom2(I3bu);
end

I3bl       = Indexset==3 & beta_l==0;
ai(I3bl)   = (x(I3bl)-lb(I3bl))./denom1(I3bl)-1;
bi(I3bl)   = phi(I3bl)./denom1(I3bl)-1;
I3bl       = Indexset==3 & beta_l~=0;
if any(I3bl)
  ai(I3bl)   = z(I3bl)./denom4(I3bl)-1;
  bi(I3bl)   = (ci(I3bl).*z(I3bl)+(di(I3bl,ones(1,n)).*DFx(I3bl,:))*z)./denom4(I3bl)-1;
end

Da(I3)     = ai(I3)+bi(I3).*ci(I3);
Db(I3)     = bi(I3).*di(I3);

I3a         = I(Indexset==3 & alpha_l==1 & alpha_u==1);
if any(I3a)
  H2(I3a,:) = bsxfun(@times,-lb(I3a)-ub(I3a)+2*x(I3a),DFx(I3a,:))+...
              2*sparse(1:length(I3a),I3a,Fx(I3a),length(I3a),n,length(I3a));
end
I3a         = I(Indexset==3 & alpha_l==1 & alpha_u~=1);
if any(I3a)
  H2(I3a,:) = bsxfun(@times,x(I3a)-lb(I3a),DFx(I3a,:))+...
              sparse(1:length(I3a),I3a,Fx(I3a),length(I3a),n,length(I3a));
end
I3a         = I(Indexset==3 & alpha_l~=1 & alpha_u==1);
if any(I3a)
  H2(I3a,:) = bsxfun(@times,x(I3a)-ub(I3a),DFx(I3a,:))+...
              sparse(1:length(I3a),I3a,Fx(I3a),length(I3a),n,length(I3a));
end

H1 = bsxfun(@times,Db,DFx);
H1 = spdiags(diag(H1)+Da,0,H1);

H  = [lambda1*H1; lambda2*H2];
