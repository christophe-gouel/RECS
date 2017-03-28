classdef recsmodelsp
% RECSMODELSP

  properties
    bounds
    dim     % Problem's dimensions {d,m,p,q}
    dima
    % FUNCTIONS - Anonymous functions that defines the model
    functions 
    infos
    nperiods
    params  % Model's parameters
    shocks
    ss      % Model's deterministic steady state
    symbols
  end

  methods
    function model = recsmodelsp(inputfiles,options)
      defaultopt = struct('Python' ,0);
      if nargin<2
        options = defaultopt;
      else
        options = catstruct(defaultopt,options);
      end

      %% Run dolo-recs on the Yaml file
      recsdirectory = fileparts(which('recsSimul'));
      inputfilesdirectory = fileparts(which(inputfiles{1}));
      recipe = fullfile(recsdirectory,'recipes','fghsubperiods.yaml');
      dolo = fullfile(recsdirectory,'Python','dolo');
      dolooptions = [' --diff --model_type=fghp --recipes="' recipe '" '];
      
      model.nperiods = length(inputfiles);
      model.functions = struct();
      model.symbols = struct();
      model.ss      = struct();
      model.ss.sss  = cell(model.nperiods,1);
      model.ss.xss  = cell(model.nperiods,1);
      
      for iperiod=1:model.nperiods
        outputfile = strrep(inputfiles{iperiod},'.yaml','model.m');

        if ispc && ~options.Python && ...
              exist(fullfile(dolo,'bin','dolo-matlab.exe'),'file')
          status = system([fullfile(dolo,'bin','dolo-matlab.exe') dolooptions ...
                           which(inputfiles{iperiod}) ' ' ...
                           fullfile(inputfilesdirectory,outputfile)]);
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
          status = system([pythonexec ' ' dolomatlab dolooptions which(inputfiles{iperiod}) ...
                           ' '  fullfile(inputfilesdirectory,outputfile)]);
        end

        if status~=0
          error('Failure to create the model file')
        end
        
        model.functions(iperiod).fname = outputfile;
      end

      modeltmp      = eval(strrep(model.functions(1).fname,'.m',''));
      model.params  = modeltmp.calibration.parameters;

      model.dim = cell(model.nperiods,4);
      for iperiod=1:model.nperiods
        modeltmp      = eval(strrep(model.functions(iperiod).fname,'.m',''));
      
        model.functions(iperiod).b  = @(s,p,varargin) OrganizeBounds(modeltmp,s,p,varargin{:});
        model.functions(iperiod).f  = modeltmp.functions.arbitrage;
        model.functions(iperiod).g  = modeltmp.functions.transition;
        model.functions(iperiod).h  = modeltmp.functions.expectation;
        
        model.infos(iperiod).incidence_matrices = modeltmp.infos.incidence_matrices;
        model.infos(iperiod).ixforward = sum(model.infos(iperiod).incidence_matrices.expectation{5},1)>=1;
   
        model.dim(iperiod,:) = {length(modeltmp.symbols.states) ...
                            length(modeltmp.symbols.controls) ...
                            length(modeltmp.symbols.expectations) ...
                            length(modeltmp.symbols.shocks)};
        
        model.ss.sss{iperiod} = modeltmp.calibration.states;
        model.ss.xss{iperiod} = modeltmp.calibration.controls;
        
        model.symbols(iperiod).states   = modeltmp.symbols.states;
        model.symbols(iperiod).controls = modeltmp.symbols.controls;
      end
        
     
      function [LB,UB,LB_s,UB_s] = OrganizeBounds(modeltmp,s,p,o)
      % ORGANIZEBOUNDS reorganizes bounds from dolo-matlab format to RECS format
        if nargin<4
          if nargout>2, o = [1; 1];
          else          o = [1; 0]; 
          end
        end
        [LB,LB_s] = modeltmp.functions.arbitrage_lb(s,p,o);
        [UB,UB_s] = modeltmp.functions.arbitrage_ub(s,p,o);        
      end
      
    end % recsmodelsp
    function SQ = StateQuant(model)
      nrep = 10000;
      nper = 200;

      SQ = cell(model.nperiods,1);
      ssim = cell(model.nperiods,1);
      for i=1:model.nperiods
        ssim{i} = zeros(nrep,model.dim{i,1},nper); 
      end
      for t=1:nper
        for i=1:model.nperiods
          if i~=1
            ssim{i}(:,:,t) = model.functions(i-1).g(ssim{i-1}(:,:,t),...
                                                    model.ss.xss{i-1}(ones(nrep,1),:),...
                                                    model.shocks{i-1}.funrand(nrep), ...
                                                    model.params);
          elseif t>1
            ssim{1}(:,:,t) = model.functions(model.nperiods).g(ssim{model.nperiods}(:,:,t-1),...
                                                              model.ss.xss{model.nperiods}(ones(nrep,1),:),...
                                                              model.shocks{model.nperiods}.funrand(nrep), ...
                                                              model.params);
          else
            ssim{1}(:,:,1) = model.ss.sss{1}(ones(nrep,1),:);
          end
        end
      end
      quant = [0 0.001 0.01 0.5 0.99 0.999 1];
      for i=1:model.nperiods
        SQ{i} = quantile(reshape(permute(real(ssim{i}(:,:,2:end)),[1 3 2]),...
                                 nrep*(nper-1),model.dim{i,1}),...
                         quant);
        if model.dim{i,1}==1, SQ{i} = SQ{i}'; end
        SQ{i} = array2table(SQ{i},...
                            'RowNames',{'min' '0.1%' '1%' 'median' '99%' '99.9%' 'max'},...
                            'VariableNames',model.symbols(i).states);
      end
    end % StateQuant
  end % methods

end % classdef