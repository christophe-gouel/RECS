classdef recsmodel
% RECSMODEL

  properties
    func    % Anonymous function that defines the model's equations
    params  % Model's parameters
    b       % Anonymous function that defines the model's bounds
    f       % Anonymous function that defines the model's equilibrium equations
    g       % Anonymous function that defines the model's transition equations
    h       % Anonymous function that defines the model's expectations definition equations
    ee      % Anonymous function that defines the model's unit-free equation errors
    e       % Values of discrete probability distribution
    w       % Probabilities of discrete probability distribution
    % FUNRAND - Random shocks generator function. Function handle that accepts an
    % integer n as input and returns a n-by-q matrix of shocks.
    funrand
    LinearSolution % Linear approximation of the rational expectations solution
    sss     % State variables at steady state
    xss     % Response variables at steady state
    zss     % Expectations variables at steady state
    dima    % Dimension of auxiliary variables
    fa      %
    ha      %
    dim     % Problem's dimensions {d,m,p}
  end % Immutable properties
  properties (Hidden=true)
    bp
    fp
    gp
    hp
    IncidenceMatrices
    ixforward % Index of response variables that appear with leads
  end % Hidden properties

  methods
    function model = recsmodel(inputfile,shocks,outputfile,options)
    % RECSMODEL Prepares a RECS model structure
    %
    % RECSMODEL uses dolo (https://github.com/albop/dolo), a Python
    % preprocessor, to convert the model described in a Yaml file to a file readable
    % by MATLAB and RECS programs. In the conversion, dolo calculates the analytic
    % representation of all partial derivatives.
    %
    % RECSMODEL(INPUTFILE) converts a model structure file, indicated by the string
    % INPUTFILE, to a m-file, readable by MATLAB and RECS programs.
    %
    % MODEL = RECSMODEL(INPUTFILE) returns MODEL a structure containing the name
    % of the model m-file and its parameters.
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
        warning('off','catstruct:DuplicatesFound')
        options = catstruct(defaultopt,options);
      end

      %% Run dolo-recs on the Yaml file
      recsdirectory = fileparts(which('recsSimul'));
      inputfiledirectory = fileparts(which(inputfile));
      dolo = fullfile(recsdirectory,'Python','dolo');

      if ispc && ~options.Python
        status = system([fullfile(dolo,'bin','dolo-recs.exe') ' ' which(inputfile) ...
                         ' ' fullfile(inputfiledirectory,outputfile)]);
      else
        dolorecs = fullfile(dolo,'bin','dolo-recs');
        setenv('PYTHONPATH',dolo)
        status = system(['python ' dolorecs ' ' which(inputfile) ...
                         ' '  fullfile(inputfiledirectory,outputfile)]);
      end

      if status~=0
        error('Failure to create the model file')
      end

      model.func        = eval(['@' strrep(outputfile,'.m','')]);
      model.params      = model.func('params');

      model.b  = @(s,p)                  model.func('b',s,[],[],[],[],[],p);
      model.f  = @(s,x,z,p,output)       model.func('f',s,x ,z ,[],[],[],p,output);
      model.g  = @(s,x,e,p,output)       model.func('g',s,x ,[],e ,[],[],p,output);
      model.h  = @(s,x,e,sn,xn,p,output) model.func('h',s,x ,[],e ,sn,xn,p,output);
      model.ee = @(s,x,z,p)              model.func('e',s,x ,z ,[],[],[],p);

      %% Incidence matrices and dimensions
      model.IncidenceMatrices  = model.func('J');
      model.dim = {size(model.IncidenceMatrices.fs,2) ...
                   size(model.IncidenceMatrices.fs,1) ...
                   size(model.IncidenceMatrices.fz,2)};
      model.ixforward = sum(model.IncidenceMatrices.hxnext,1)>=1;


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

        %% Find steady state
        [sss0,xss0] = model.func('ss');
        if ~isempty(sss0) && ~isempty(xss0)
          [sss,xss,zss,exitflag] = recsSS(model,sss0,xss0,options);
          if exitflag==1
            model.sss = sss;
            model.xss = xss;
            model.zss = zss;
          end
        end
      end % shocks and steady state
    end % recsmodel
    function SQ = StateQuant(model)
      nrep = 10000;
      nper = 200;
      d = model.dim{1};
      ssim = zeros(nrep,d,nper+1);
      ssim(:,:,1) = model.sss(ones(nrep,1),:);
      output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
      for t=1:nper
        ssim(:,:,t+1) = model.g(ssim(:,:,t),model.xss(ones(nrep,1),:),...
                                model.funrand(nrep),model.params,output);
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