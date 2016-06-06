function [interp,X,Z,exitflag,output] = recsFirstGuessSP(model,interp,options)
% RECSFIRSTGUESSSP finds a first guess using the perfect foresight solution
%
% RECSFIRSTGUESSSP tries to find a first-guess for the rational expectations model
% by solving the corresponding perfect foresight solution on all the grid points
% of the state variables. By default, it considers that the model goes back to
% its deterministic steady state in 50 periods.
%
% INTERP = RECSFIRSTGUESSSP(MODEL,INTERP) uses the interpolation structure
% defined in the structure INTERP to fit the perfect foresight solution of the
% model defined in the structure MODEL. RECSFIRSTGUESSSP returns an
% interpolation structure, INTERP, containing the first guess.
% INTERP is a structure, which has to include the following fields:
%    fspace       : a definition structure for the interpolation family (created
%                   by the function fundef)
% MODEL is an object created by recsmodel.
%
% INTERP = RECSFIRSTGUESSSP(MODEL,INTERP,OPTIONS) solves the problem
% with the parameters defined by the structure OPTIONS. The fields of the
% structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver (default:
%                       empty structure)
%    fgmethod         : 'auto' (default), 'perturbation', 'perfect-foresight' or
%                       'steady-state'
%    T                : integer defining the time horizon at which the model is
%                       supposed to converge to its steady state (default: 50)
%
% [INTERP,X] = RECSFIRSTGUESSSP(MODEL,INTERP) returns X, n-by-m matrix,
% containing the value of the response variables in the first period.
%
% [INTERP,X,Z] = RECSFIRSTGUESSSP(MODEL,INTERP) returns Z, n-by-p matrix,
% containing the value of the expectations variables in the first period.
%
% [INTERP,X,Z,EXITFLAG] = RECSFIRSTGUESSSP(MODEL,INTERP)
%
% [INTERP,X,Z,EXITFLAG,OUTPUT] = RECSFIRSTGUESSSP(MODEL,INTERP)
%
% See also RECSSOLVEDETERMINISTICPBSP, RECSSOLVEREESP, RECSSSSP.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct('fgmethod','auto',...
                    'T'       ,50);
if nargin <=2
  options = defaultopt;
else
  options = catstruct(defaultopt,options);
end

s        = interp.s;

dim      = model.dim;
nperiods = model.nperiods;

fgmethod = options.fgmethod;
T        = options.T;

%% Solve for the deterministic steady state
[sss,xss,zss] = recsSSSP(model,model.ss.sss,model.ss.xss,...
                         catstruct(options,struct('display',0)));

%% Find first-guess
if strcmpi(fgmethod,'auto'), fgmethod = 'perfect-foresight'; end

switch fgmethod
  case 'perfect-foresight'
    %% Solve the perfect foresight problem on each point of the grid
    X        = cellfun(@(S,dimX) zeros(size(S,1),dimX),s,dim(:,2),'UniformOutput', false);
    Z        = cellfun(@(S,dimX) zeros(size(S,1),dimX),s,dim(:,3),'UniformOutput', false);
    exitflag = cellfun(@(S) zeros(size(S,1),1),s,'UniformOutput', false);
    N        = cellfun(@(S) zeros(size(S,1),1),s,'UniformOutput', false);

    parfor i=1:nperiods
      for j=1:size(s{i},1)
        [x,~,z,~,exitflag{i}(j),N{i}(j)] = recsSolveDeterministicPbSP(model,...
                                                          s{i}(j,:),i,...
                                                          T,xss,zss,sss,options);
        X{i}(j,:) = x{i}(1,:);
        Z{i}(j,:) = z{i}(1,:);
      end
    end
    output = struct('exitflag',exitflag,...
                    'N'       ,N);
    exitflag = ~any(cat(1,exitflag{:})~=1);

  case 'steady-state'
    %% Response variables equal to their steady-state values
    LB = cell(nperiods,1);
    UB = cell(nperiods,1);
    for i=1:nperiods
      [LB{i},UB{i}] = model.functions(i).b(s{i},model.params);
    end
    X = cellfun(@(XSS,S,LB,UB) min(max(repmat(XSS,size(S,1),1),LB),UB),...
                xss,s,LB,UB);
    Z = cellfun(@(S,dimZ) NaN(size(S,1),dimZ),s,dim(:,3));
    exitflag = 1;
    output   = [];

end % switch fgmethod

%% Prepare output
interp.cX = cellfun(@funfitxy,interp.fspace,interp.Phi,X,'UniformOutput', false);
interp.cZ = cellfun(@funfitxy,interp.fspace,interp.Phi,Z,'UniformOutput', false);
interp.X  = X;
interp.Z  = Z;

