classdef recsmodelsp
% RECSMODELSP

  properties
    params  % Model's parameters
    % FUNCTIONS - Anonymous functions that defines the model
    functions 
    dim     % Problem's dimensions {d,m,p,q}
    nperiods
    shocks
    bounds
    infos
  end

  methods
    function model = recsmodelsp(inputfiles,options)
      defaultopt = struct('Python' ,0);
      if nargin<2
        options = defaultopt;
      else
        warning('off','catstruct:DuplicatesFound')
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
      
      for iperiod=1:model.nperiods
        outputfile = strrep(inputfiles{iperiod},'.yaml','model.m');

        if ispc && ~options.Python && ...
              exist(fullfile(dolo,'bin','dolo-matlab.exe'),'file')
          status = system([fullfile(dolo,'bin','dolo-matlab.exe') dolooptions ...
                           which(inputfiles{iperiod}) ' ' ...
                           fullfile(inputfilesdirectory,outputfile)]);
        else
          dolomatlab = fullfile(dolo,'bin','dolo-matlab');
          setenv('PYTHONPATH',dolo)
          status = system(['python ' dolomatlab dolooptions which(inputfiles{iperiod}) ...
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
  end % methods

end % classdef