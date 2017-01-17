classdef recsmodel
% RECSMODEL

  properties
    dim     % Problem's dimensions {d,m,p,q}
    dima    % Dimension of auxiliary variables
    func    % Anonymous function that defines the model's equations
    params  % Model's parameters
    % FUNCTIONS - Anonymous functions that defines the model
    functions
    % FUNRAND - Random shocks generator function. Function handle that accepts an
    % integer n as input and returns a n-by-q matrix of shocks.
    funrand
    LinearSolution % Linear approximation of the rational expectations solution
    sss     % State variables at steady state
    xss     % Response variables at steady state
    zss     % Expectations variables at steady state
    shocks
    symbols % Symbols of parameters, shocks, and variables
    infos   % Various information about the model
  end

  methods
    function model = recsmodel(inputfile,shocks,outputfile,options)
    % RECSMODEL Prepares a recsmodel object
    %
    % RECSMODEL uses dolo (https://github.com/albop/dolo), a Python
    % preprocessor, to convert the model described in a Yaml file to a file readable
    % by MATLAB and RECS programs. In the conversion, dolo calculates the analytic
    % representation of all partial derivatives.
    %
    % MODEL = RECSMODEL(INPUTFILE) converts a model Yaml file, indicated by the
    % string INPUTFILE, to a m-file, readable by MATLAB and RECS programs. It
    % returns MODEL an object containing the name of the model m-file, its
    % parameters, and other properties.
    %
    % MODEL = RECSMODEL(INPUTFILE,SHOCKS) prepares in MODEL the shocks
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
    % MODEL = RECSMODEL(INPUTFILE,SHOCKS,OUTPUTFILE) uses the string OUTPUTFILE
    % to name the m-file containing the model.
    %
    % MODEL = RECSMODEL(INPUTFILE,SHOCKS,OUTPUTFILE,OPTIONS) uses the options
    % defined by the structure OPTIONS. The fields of the structure are
    %    display          : 1 to display the steady state if found (default: 1)
    %    eqsolver         : solver used to find the steady state 'fsolve', 'lmmcp'
    %                       (default), 'ncpsolve' or 'path'
    %    eqsolveroptions  : options structure to be passed to eqsolver
    %    Python           : 1 to call Python directly instead of the executable file
    %                       (default: 0, only for Windows and for developement)
      if nargin<3 || isempty(outputfile)
        outputfile = strrep(inputfile,'.yaml','model.m');
      end

      defaultopt = struct('Python' ,0);
      if nargin<4
        options = defaultopt;
      else
        options = catstruct(defaultopt,options);
      end

      %% Run dolo-recs on the Yaml file
      recsdirectory = fileparts(which('recsSimul'));
      inputfiledirectory = fileparts(which(inputfile));
      dolo = fullfile(recsdirectory,'Python','dolo');
      dolooptions = ' --diff  --model_type=fgh2 ';

      if ispc && ~options.Python && ...
            exist(fullfile(dolo,'bin','dolo-matlab.exe'),'file')
        status = system([fullfile(dolo,'bin','dolo-matlab.exe') dolooptions ...
                         which(inputfile) ' ' ...
                         fullfile(inputfiledirectory,outputfile)]);
      else
        PythonVirtualEnvdirectory = fullfile(recsdirectory,'Python','PythonVirtualEnv');
        if ispc
          if exist(fullfile(PythonVirtualEnvdirectory,'python.exe'),'file')
            pythonexec = fullfile(PythonVirtualEnvdirectory,'python.exe');
          elseif exist(fullfile(PythonVirtualEnvdirectory,'bin', ...
                                'python.exe'),'file')
            pythonexec = fullfile(PythonVirtualEnvdirectory,'bin', ...
                                  'python.exe');
          elseif exist(fullfile(PythonVirtualEnvdirectory,'Scripts', ...
                                'python.exe'),'file')
            pythonexec = fullfile(PythonVirtualEnvdirectory,'Scripts', ...
                                  'python.exe');
          else
            pythonexec = 'python.exe';
          end
        else
          if exist(fullfile(PythonVirtualEnvdirectory,'bin','python'),'file')
            pythonexec = fullfile(PythonVirtualEnvdirectory,'bin','python');
          else
            pythonexec = 'python';
          end
        end % ispc
        dolomatlab = fullfile(dolo,'bin','dolo-matlab');
        setenv('PYTHONPATH',dolo)
        if ~exist(pythonexec,'file'), pythonexec = 'python'; end
        status = system([pythonexec ' ' dolomatlab dolooptions which(inputfile) ...
                         ' '  fullfile(inputfiledirectory,outputfile)]);
      end

      if status~=0
        error('Failure to create the model file')
      end

      modeltmp      = eval(strrep(outputfile,'.m',''));
      model.params  = modeltmp.calibration.parameters;
      model.symbols = modeltmp.symbols;

      sss0 = modeltmp.calibration.states;
      xss0 = modeltmp.calibration.controls;

      model.functions.b  = @(s,p,varargin) OrganizeBounds(modeltmp,s,p,varargin{:});
      model.functions.f  = modeltmp.functions.arbitrage;
      model.functions.g  = modeltmp.functions.transition;
      model.functions.h  = modeltmp.functions.expectation;
      model.functions.ee = @(s,x,z,p) NaN;

      %% Incidence matrices and dimensions
      IM = modeltmp.infos.incidence_matrices;
      model.infos.IncidenceMatrices = struct('fs',IM.arbitrage{1},...
                                             'fx',IM.arbitrage{2},...
                                             'fz',IM.arbitrage{3},...
                                             'gs',IM.transition{1},...
                                             'gx',IM.transition{2},...
                                             'ge',IM.transition{3},...
                                             'hs',IM.expectation{1},...
                                             'hx',IM.expectation{2},...
                                             'he',IM.expectation{3},...
                                             'hsnext',IM.expectation{4},...
                                             'hxnext',IM.expectation{5},...
                                             'lbs',IM.arbitrage_lb{1},...
                                             'ubs',IM.arbitrage_ub{1});

      model.dim = {size(model.infos.IncidenceMatrices.fs,2) ...
                   size(model.infos.IncidenceMatrices.fs,1) ...
                   size(model.infos.IncidenceMatrices.fz,2) ...
                   length(model.symbols.shocks)};
      model.infos.ixforward = sum(model.infos.IncidenceMatrices.hxnext,1)>=1; % Index of response variables that appear with leads

      %% Identify variables for which bounds are variable
      model.infos.ixvarbounds = logical([sum(model.infos.IncidenceMatrices.lbs,2)...
                          sum(model.infos.IncidenceMatrices.ubs,2)]);
      model.infos.nxvarbounds = int16(sum(model.infos.ixvarbounds,1));
      if any(model.infos.nxvarbounds)
        model    = mcptransform(model);
      else
        model.functions.bp = @(s,p,varargin) OrganizeBounds(modeltmp,s,p,varargin{:});;
        model.functions.fp = @(s,x,w,v,z,p,o) modeltmp.functions.arbitrage(s,x,z,p,o);
        model.functions.gp = modeltmp.functions.transition;
        model.functions.hp = modeltmp.functions.expectation;
      end

      %% Identify model type
      if all(~[model.infos.IncidenceMatrices.hs(:); ...
               model.infos.IncidenceMatrices.hx(:); ...
               model.infos.IncidenceMatrices.he(:)])
        model.infos.model_type = 'fgh1';
      else
        model.infos.model_type = 'fgh2';
      end

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
        [model.shocks.e,model.shocks.w] = qnwnorm(order,Mu,Sigma);
        % Random number generator
        model.shocks.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,q)*R;

        %% Find steady state
        if ~isempty(sss0) && ~isempty(xss0)
          [sss,xss,zss,exitflag] = recsSS(model,sss0,xss0,options);
          if exitflag==1
            model.sss = sss;
            model.xss = xss;
            model.zss = zss;
          end
        end
      end % shocks and steady state

      %% Equation type
      model.infos.eq_type = 'mcp';
      if all(~[model.infos.IncidenceMatrices.lbs(:); ...
               model.infos.IncidenceMatrices.ubs(:)])
        if ~isempty(sss0)
          [LB,UB] = model.functions.b(sss0,model.params);
          if all(isinf([LB(:); UB(:);]))
            model.infos.eq_type = 'cns';
          end
        end
      end

      function [LB,UB,LB_s,UB_s] = OrganizeBounds(modeltmp,s,p,o)
      % ORGANIZEBOUNDS reorganizes bounds from dolo-matlab format to RECS format
        if nargin<4
          if nargout>2, o = [1; 1];
          else,         o = [1; 0];
          end
        end
        [LB,LB_s] = modeltmp.functions.arbitrage_lb(s,p,o);
        [UB,UB_s] = modeltmp.functions.arbitrage_ub(s,p,o);
      end

    end % recsmodel
    function SQ = StateQuant(model)
      nrep = 10000;
      nper = 200;
      d = model.dim{1};
      ssim = zeros(nrep,d,nper+1);
      ssim(:,:,1) = model.sss(ones(nrep,1),:);
      for t=1:nper
        ssim(:,:,t+1) = model.functions.g(ssim(:,:,t),...
                                          model.xss(ones(nrep,1),:),...
                                          model.shocks.funrand(nrep),model.params);
      end
      quant = [0 0.001 0.01 0.5 0.99 0.999 1];
      SQ = quantile(reshape(permute(real(ssim(:,:,2:end)),[1 3 2]),nrep*nper,d),...
                     quant);
      if d==1, SQ = [quant' SQ'];
      else     SQ = [quant' SQ ];
      end
    end % StateQuant
  end % methods

end % classdef