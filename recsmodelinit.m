function model = recsmodelinit(inputfile,shocks,outputfile)
% RECSMODELINIT Prepares a RECS model structure
%
% RECSMODELINIT uses dolo (https://github.com/albop/dynare-python), a python
% preprocessor, to convert the model described in a yaml file to a file readable by
% MATLAB and RECS programs. In the conversion, dolo calculates the analytic
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

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin<3, outputfile = strrep(inputfile,'.yaml','model.m'); end

%% Run dolo-recs on the yaml file
recsdirectory = fileparts(which('recsSimul'));
inputfiledirectory = fileparts(which(inputfile));

if ~(ispc || strcmp(computer('arch'),'glnx86'))
  error('Not available on this platform')
end

if ispc
  txt = 'win32';
  extension = '.exe';
else
  txt = computer('arch');
  extension = '';
end
status = system([fullfile(recsdirectory,'exe',txt,['dolo-recs' extension]) ...
                 ' ' which(inputfile) ' ' fullfile(inputfiledirectory,outputfile)]);

if status~=0
  error('Failure to create the model file')
end

%% Pack model structure
model.func        = eval(['@' strrep(outputfile,'.m','')]);
model.params      = model.func('params');

% Prepare shocks information
if nargin>=2 && ~isempty(shocks)
  % Unpack shocks
  Mu    = shocks.Mu;
  q     = length(Mu);
  order = shocks.order;
  if isscalar(order) && q>1, order = order*ones(1,q); end
  Sigma = shocks.Sigma;

  % Check if Sigma is positive-semidefinite and symmetric
  if ~isequal(Sigma,Sigma') || ...
        ~all(diag(Sigma)>=0) || ...
        ~all(diag(Sigma)'>=sum(abs(Sigma-diag(diag(Sigma)))))
    error('shocks.Sigma is not a symmetric, positive-semidefinite matrix')
  end
  R = chol(Sigma);

  % Gaussian quadrature
  [model.e,model.w] = qnwnorm(order,Mu,Sigma);
  % Random number generator
  model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,q)*R;
end