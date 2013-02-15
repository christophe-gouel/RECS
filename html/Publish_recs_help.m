function Publish_recs_help(website)
% PUBLISH_RECS_HELP publishes help pages to html

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin < 1, website = 0; end

recsdirectory   = fileparts(which('recsSimul'));
if website
  targetdirectory = fullfile(recsdirectory,'www');
else
  targetdirectory = fullfile(recsdirectory,'html');
end
PublishOptions = struct('outputDir',targetdirectory);
if exist('html.xsl','file')
  PublishOptions = catstruct(PublishOptions,struct('stylesheet','html.xsl'));
  if website
    PublishOptions = catstruct(PublishOptions,struct('stylesheet','htmlwithmenu.xsl'));
  end
end
PublishOptionsNoExec = catstruct(PublishOptions,struct('evalCode',false));
PublishOptionsNoShow = catstruct(PublishOptions,struct('showCode',false));

%% Clear the target directory
delete(fullfile(targetdirectory,'*.png'));
delete(fullfile(targetdirectory,'*.txt'));
delete(fullfile(targetdirectory,'*.yaml'));
delete(fullfile(targetdirectory,'*.html'));

%% Documentation
publish('recs_product_page.m',PublishOptions);
if website
  copyfile(fullfile(targetdirectory,'recs_product_page.html'),...
           fullfile(targetdirectory,'index.html'));
end

% Getting started
publish('getting_started.m',PublishOptions);
publish('installation.m',PublishOptions);
publish('tutorial.m',PublishOptions);
publish('def_sre.m',PublishOptions);
publish('MCP.m',PublishOptions);

% User guide
publish('user_guide.m',PublishOptions);
publish('ug_setting_up.m',PublishOptions);
publish('ug_model_files.m',PublishOptions);
addpath(fullfile(recsdirectory,'demos'))
publish('ug_model_struct.m',PublishOptions);
publish('ug_interpolation.m',PublishOptions);
publish('ss.m',PublishOptions);
publish('first_guess.m',PublishOptions);
publish('solve_REE.m',PublishOptions);
publish('simulate.m',PublishOptionsNoShow);
rmpath(fullfile(recsdirectory,'demos'))
publish('calibration.m',PublishOptions);
publish('accuracy.m',PublishOptions);
publish('finite_horizon.m',PublishOptions);
publish('deterministic.m',PublishOptions);
publish('ug_solvers_eq.m',PublishOptions);
publish('ug_methods.m',PublishOptions);
publish('ug_options.m',PublishOptions);

% Others
publish('recs_functions.m',PublishOptions);
publish('demos.m',PublishOptions);
publish('pathnotinstalled.m',PublishOptions);

%% Functions
if website
  FunctionList = {'recsAccuracy',...
                  'recsAuxiliary',...
                  'recsCheck',...
                  'recsConvert',...
                  'recsFirstGuess',...
                  'recsinterpinit',...
                  'recsmodelinit',...
                  'recsSimul',...
                  'recsSolveDeterministicPb',...
                  'recsSolveREE',...
                  'recsSolveREEFiniteHorizon',...
                  'recsSS'};
  for fn = FunctionList
    fid = fopen(fullfile(targetdirectory,[fn{1} '.html']),'w');
    fprintf(fid,'%s',help2html([fn{1} '.m']));
    fclose(fid);
  end
end

%% License
copyfile(fullfile(recsdirectory,'LICENSE.txt'),fullfile(targetdirectory,'LICENSE.txt'));

%% Demonstration
currentfolder = cd(fullfile(recsdirectory,'demos'));
copyfile('cs1.yaml',fullfile(targetdirectory,'cs1.txt'));
copyfile('gro1.yaml',fullfile(targetdirectory,'gro1.txt'));
copyfile('gro2.yaml',fullfile(targetdirectory,'gro2.txt'));
copyfile('sto1.yaml',fullfile(targetdirectory,'sto1.txt'));
copyfile('sto2.yaml',fullfile(targetdirectory,'sto2.txt'));
copyfile('sto4.yaml',fullfile(targetdirectory,'sto4.txt'));
copyfile('sto5.yaml',fullfile(targetdirectory,'sto5.txt'));
copyfile('sto6.yaml',fullfile(targetdirectory,'sto6.txt'));
publish('clearpublish.m',PublishOptions);
publish('cs1.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('cs2.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('gro1.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('gro2.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto1.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto2.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto3.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto4.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto5.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto6.m',PublishOptions);
publish('cs1model.m',PublishOptionsNoExec);
publish('gro1model.m',PublishOptionsNoExec);
publish('gro2model.m',PublishOptionsNoExec);
publish('sto1model.m',PublishOptionsNoExec);
publish('sto2model.m',PublishOptionsNoExec);
publish('sto4model.m',PublishOptionsNoExec);
publish('sto5model.m',PublishOptionsNoExec);
publish('sto6model.m',PublishOptionsNoExec);
delete(fullfile(targetdirectory,'clearpublish.html'));
cd(currentfolder)

%% Build search database
if ~website, builddocsearchdb(targetdirectory); end