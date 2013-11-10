function [x,s,z,F,exitflag,N] = recsSolveDeterministicPb(model,s0,T,xss,zss,sss,options)
% RECSSOLVEDETERMINISTICPB Solves a perfect foresight problem
%
% X = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS) tries to find the perfect
% foresight solution of the model defined in the object MODEL. The initial
% values of state variable are provided by the 1-by-d vector S0. Time horizon
% (number of periods) is given by the integer T. XSS, ZSS and SSS are,
% respectively, a 1-by-m vector containing the values of response variables at the
% deterministic steady state, a 1-by-p vector containing the values of expectations
% variables at steady state, and a 1-by-d vector containing the values of the state
% variables at steady state. RECSSOLVEDETERMINISTICPB returns X, a T-by-m matrix,
% the value of response variables over the time horizon.
%
% X = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,OPTIONS) solves the problem with the
% parameters defined by the structure OPTIONS. The fields of the structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver (default:
%                       empty structure)
%
% [X,S] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns S, a T-by-d
% matrix, containing the value of state variables over the time horizon.
%
% [X,S,Z] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns Z, a
% T-by-p matrix, containing the value of expectations variables over the time
% horizon.
%
% [X,S,Z,F] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns F, a
% T-by-m matrix, containing the values of equilibrium equations over the time
% horizon.
%
% [X,S,Z,F,EXITFLAG] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...)
%
% [X,S,Z,F,EXITFLAG,N] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...)
%
% See also RECSFIRSTGUESS, RECSSOLVEREE, RECSSS, SCP.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                      ...
    'checkfinalstate' , 0                                ,...
    'eqsolver'        , 'lmmcp'                          ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off' ,...
                               'DerivativeCheck', 'off' ,...
                               'Jacobian'       , 'on'));
if nargin <=6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;

[d,m,p] = model.dim{1:3};
n       = T*(m+p+d);

e      = model.shocks.w'*model.shocks.e;
fp     = model.functions.fp;
gp     = model.functions.gp;
hp     = model.functions.hp;
params = model.params;

nx = model.infos.nxvarbounds;
M  = m+sum(nx);

n = n+T*sum(nx);

w = zeros(1,nx(1));
v = zeros(1,nx(2));

[LBx,UBx] = model.functions.bp(sss,params);
LB = [LBx -inf(1,p+d)];
UB = [UBx +inf(1,p+d)];

LB = reshape(LB(ones(T,1),:)',n,1);
UB = reshape(UB(ones(T,1),:)',n,1);

X = [xss w v zss sss];
X = X(ones(T,1),:)';
X = reshape(X,n,1);

ix = 1:M;
iz = (M+1):(M+p);
is = (M+p+1):(M+p+d);

vec     = @(X) X(:);
iX2iXT  = @(iX,Tperiod,dimX) vec((repmat(iX,length(Tperiod),1)+...
                                  (M+p+d)*repmat(Tperiod',1,dimX))');
izprevT = iX2iXT(iz,0:T-2,p);
isprevT = iX2iXT(is,0:T-2,d);
ixT     = iX2iXT(ix,0:T-1,M);
izT     = iX2iXT(iz,0:T-1,p);
isT     = iX2iXT(is,0:T-1,d);
ixnextT = iX2iXT(ix,1:T-1,M);
iznextT = iX2iXT(iz,1:T-1,p);
isnextT = iX2iXT(is,1:T-1,d);

%% Number of nonzero elements in the Jacobian
IM = structfun(@(X) sum(X(:)),...
               model.infos.IncidenceMatrices,'UniformOutput',false);
nnzJac = T*(IM.fx+2*(IM.lbs+IM.ubs)+IM.fz+IM.hx+p+IM.hsnext+IM.gx+d)+... % Main diagonal blocks
         (T-1)*(IM.fs+IM.lbs+IM.ubs+IM.hs+IM.gs)+...                     %   Subdiagonal blocks
         (T-1)*IM.hxnext;                                                % Superdiagonal blocks
eqsolveroptions.nnzJ = nnzJac;

%% Solve deterministic problem

% Precalculation of indexes for the sparse Jacobian
tmp      = zeros(M+p+d,M+p+d,T);
[~,grid] = blktridiag(tmp,tmp(:,:,1:end-1),tmp(:,:,1:end-1));

SCPSubProblem = @(X0,S0) runeqsolver(@recsDeterministicPb,X0,LB,UB,eqsolver, ...
                                     eqsolveroptions,fp,gp,hp,S0,xss,p,e, ...
                                     params,nx,grid,nnzJac,izprevT,isprevT,ixT,izT, ...
                                     isT,ixnextT,iznextT,isnextT);

% Simple continuation problem applied on a Newton solve
[X,F,exitflag,N] = SCP(X,s0,sss,SCPSubProblem,1);
if exitflag~=1
  warning('RECS:FailureDeterministic',...
          'Failure to find the perfect foresight solution');
end

%% Prepare output
X = reshape(X,M+p+d,T)';
x = X(:,1:m);
z = X(:,(M+1):(M+p));
s = [s0; X(1:end-1,(M+p+1):(M+p+d))];
if ~isempty(F), F = reshape(F,M+p+d,T)'; end

%% Check is the final state is self-replicating (meaning that sT=sss and xT=xss)
if options.checkfinalstate
  deltafinal = max(abs([s(end,:)-sss x(end,:)-xss]));
  if deltafinal>sqrt(eps)
    warning('RECS:FinalState',...
            'Final state is not self-replicating (max(|delta|)=%g)',...
            deltafinal);
  end
end

