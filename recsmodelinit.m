function [model,statequant] = recsmodelinit(inputfile,shocks,outputfile,options)
% RECSMODELINIT Prepares a RECS model structure
%
% RECSMODELINIT uses dolo (https://github.com/albop/dolo), a Python
% preprocessor, to convert the model described in a Yaml file to a file readable
% by MATLAB and RECS programs. In the conversion, dolo calculates the analytic
% representation of all partial derivatives.
%
% RECSMODELINIT(INPUTFILE) converts a model structure file, indicated by the string
% INPUTFILE, to a m-file, readable by MATLAB and RECS programs.
%
% MODEL = RECSMODELINIT(INPUTFILE) returns MODEL a structure containing the name
% of the model m-file and its parameters.
%
% MODEL = RECSMODELINIT(INPUTFILE,SHOCKS) prepares in MODEL the shocks
% information by using the structure SHOCKS. The fields of the SHOCKS define the
% parameters of a multivariate normal distribution and its approximation by
% gaussian quadrature. They are
%    Mu    : 1-by-q vector of mean
%    order : 1-by-q vector (or scalar) defining the number of nodes of each
%            shock variables in the gaussian quadrature. If a scalar is passed,
%            it is extended to allow the same number of nodes for all variables.
%    Sigma : q-by-q, symmetric, positive-semidefinite, covariance matrix
% If the shocks do not follow a multivariate normal distribution, the shocks
% information has to be produce manually.
%
% MODEL = RECSMODELINIT(INPUTFILE,SHOCKS,OUTPUTFILE) uses the string OUTPUTFILE
% to name the m-file containing the model.
%
% MODEL = RECSMODELINIT(INPUTFILE,SHOCKS,OUTPUTFILE,OPTIONS) uses the options
% defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to display the steady state if found (default: 1)
%    eqsolver         : solver used to find the steady state 'fsolve', 'lmmcp'
%                       (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [MODEL,STATEQUANT] = RECSMODELINIT(INPUTFILE,SHOCKS,...) returns the
% 6-by-d matrix STATEQUANT that contains the quantiles (0%, 0.1%, 1%,
% 99%, 99.9% 100%) for the state variables when maintaining response
% variables at their values at the deterministic steady
% state. STATEQUANT is useful for determining the extent of the
% state space for exogenous state variables.
%
% See also RECSINTERPINIT.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<3 || isempty(outputfile)
  outputfile = strrep(inputfile,'.yaml','model.m');
end

defaultopt = struct(                                      ...
    'display'           , 1                              ,...
    'eqsolver'          , 'lmmcp'                        ,...
    'eqsolveroptions'   , struct('DerivativeCheck','off'));
if nargin<4
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

%% Run dolo-recs on the Yaml file
recsdirectory = fileparts(which('recsSimul'));
inputfiledirectory = fileparts(which(inputfile));
dolo = fullfile(recsdirectory,'Python','dolo','src');

if ispc
  status = system([fullfile(dolo,'bin','dolo-recs.exe') ' ' which(inputfile) ...
                   ' ' fullfile(inputfiledirectory,outputfile)]);
elseif isunix && ~ismac
  dolorecs = fullfile(dolo,'bin','dolo-recs.py');
  status = system(['PYTHONPATH=' dolo ' python ' dolorecs ' ' which(inputfile) ...
                   ' '  fullfile(inputfiledirectory,outputfile)]);
else
  error('Not available on this platform')
end

if status~=0
  error('Failure to create the model file')
end

%% Pack model structure
model.func        = eval(['@' strrep(outputfile,'.m','')]);
model.params      = model.func('params');

%% Prepare shocks information & find steady state
if nargin>=2 && ~isempty(shocks)
  % Unpack shocks
  Mu    = shocks.Mu;
  q     = length(Mu);
  order = shocks.order;
  if isscalar(order) && q>1, order = order*ones(1,q); end
  Sigma = shocks.Sigma;
  if length(Sigma)==length(Sigma(:)), Sigma = diag(Sigma); end

  % Check if Sigma is positive-semidefinite and symmetric
  if ~isequal(Sigma,Sigma') || ...
        ~all(diag(Sigma)>=0) || ...
        ~all(diag(Sigma)'>=sum(abs(Sigma-diag(diag(Sigma)))))
    error('shocks.Sigma is not a symmetric, positive-semidefinite matrix')
  end
  if isequal(diag(diag(Sigma)),Sigma)
    % Diagonal matrix
    R = sqrt(Sigma);
  else
    R = chol(Sigma);
  end

  % Gaussian quadrature
  [model.e,model.w] = qnwnorm(order,Mu,Sigma);
  % Random number generator
  model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,q)*R;

  %% Check steady state
  [sss0,xss0] = model.func('ss');
  if ~isempty(sss0) && ~isempty(xss0)
    [sss,xss,zss,exitflag] = recsSS(model,sss0,xss0,options);
    if exitflag==1 
      model.sss = sss;
      model.xss = xss;
      model.zss = zss;
      if options.display==1
        disp('Deterministic steady state')
        fprintf(1,' State variables:\n\t\t')
        fprintf(1,'%0.4g\t',sss)
        fprintf(1,'\n\n Response variables:\n\t\t')
        fprintf(1,'%0.4g\t',xss)
        fprintf(1,'\n\n Expectations variables:\n\t\t')
        fprintf(1,'%0.4g\t',zss)
        fprintf(1,'\n\n')
      end
    end
  end

  %% Simulate state variables dynamics with response at steady-state
  % Provide a first guess of the relevant state space
  if nargout==2 && exitflag==1
    nrep        = 10000;
    nper        = 200;
    d           = size(sss,2);
    ssim        = zeros(nrep,d,nper+1);
    ssim(:,:,1) = sss(ones(nrep,1),:);
    output      = struct('F',1,'Js',0,'Jx',0,'Jz',0);
    for t=1:nper
      ssim(:,:,t+1) = model.func('g',ssim(:,:,t),xss(ones(nrep,1),:),[],...
                                 model.funrand(nrep),[],[],model.params,output);
    end
    statequant = quantile(reshape(permute(ssim(:,:,2:end),[1 3 2]),nrep*nper,d),...
                          [0 0.001 0.01 0.99 0.999 1]);
  end
end