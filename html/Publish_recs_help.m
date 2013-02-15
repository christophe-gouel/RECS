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
if website
  delete(fullfile(targetdirectory,'*.m'));
end

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
  txt = fileread(fullfile(targetdirectory,'recs_functions.html'));
  % Replace matlab:doc
  pattern = '"matlab:doc\(''(\w*)''\)"';
  txt = regexprep(txt,pattern,'"$1.html"');
  fid = fopen(fullfile(targetdirectory,'recs_functions.html'),'w');
  fprintf(fid,'%s',txt);
  fclose(fid);
  
  FunctionList = {'recsAccuracy',...
                  'recsAuxiliary',...
                  'recsCheck',...
                  'recsConvert',...
                  'recsdemos',...
                  'recsFirstGuess',...
                  'recsinterpinit',...
                  'recsmodelinit',...
                  'recsSimul',...
                  'recsSolveDeterministicPb',...
                  'recsSolveREE',...
                  'recsSolveREEFiniteHorizon',...
                  'recsSS',...
                  'lmmcp',...
                  'nsoli',...,
                  'recspathmcp',...
                  'SA'};
  for fn = FunctionList
    copyfile(fullfile(recsdirectory,[fn{1} '.m']),fullfile(targetdirectory,[fn{1} '.m']));
    txt = help2html([fn{1} '.m']);
    % Replace see also
    pattern = '"matlab:helpwin (\w*)"';
    txt = regexprep(txt,pattern,'"$1.html"');
    % Replace view code
    pattern = '"matlab:edit (\w*).m"';
    txt = regexprep(txt,pattern,'"$1.m"');
    % Replace css file
    pattern = '"file.*\.css"';
    txt = regexprep(txt,pattern,'"layoutfunctionhelp.css"');
    % Replace Matlab File help
    pattern = 'MATLAB File Help';
    txt = regexprep(txt,pattern,'Function Help');
    % Replace Default Topics
    pattern = '"matlab:helpwin">Default Topics';
    txt = regexprep(txt,pattern,'"recs_functions.html">Function Reference');
    fid = fopen(fullfile(targetdirectory,[fn{1} '.html']),'w');
    fprintf(fid,'%s',txt);
    fclose(fid);
  end
end

%% License
copyfile(fullfile(recsdirectory,'LICENSE.txt'),fullfile(targetdirectory,'LICENSE.txt'));

%% Demonstration
currentfolder = cd(fullfile(recsdirectory,'demos'));
DemoFileList = {'cs1','cs2','gro1','gro2','sto1','sto2','sto3','sto4', ...
                'sto5','sto6'};
for demo=DemoFileList
  publish('clearpublish.m',PublishOptions);
  publish([demo{1} '.m'],PublishOptions);
end
YamlFileList = {'cs1','gro1','gro2','sto1','sto2','sto4','sto5','sto6'};
for yaml=YamlFileList
  copyfile([yaml{1} '.yaml'],fullfile(targetdirectory,[yaml{1} '.txt']));
  publish([yaml{1} 'model.m'],PublishOptionsNoExec);
end
delete(fullfile(targetdirectory,'clearpublish.html'));
cd(currentfolder)

%% Build search database
if ~website, builddocsearchdb(targetdirectory); end