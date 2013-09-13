function Publish_recs_help(website)
% PUBLISH_RECS_HELP publishes help pages to html

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin < 1, website = 0; end

spmd, warning('off','backtrace'); end

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
if website, delete(fullfile(targetdirectory,'*.m')); end

%% Documentation
publish('recs_product_page.m',PublishOptions);
if website
  copyfile(fullfile(targetdirectory,'recs_product_page.html'),...
           fullfile(targetdirectory,'index.html'));
end

HelpFileList = {'getting_started.m',...
                'installation.m',...
                'tutorial.m',...
                'def_sre.m',...
                'MCP.m',...
                'solution_methods.m',...
                'user_guide.m',...
                'ug_setting_up.m',...
                'ug_model_files.m',...
                'ug_model_struct.m',...
                'ug_interpolation.m',...
                'ss.m',...
                'first_guess.m',...
                'solve_REE.m',...
                'simulate.m',...
                'accuracy.m',...
                'finite_horizon.m',...
                'deterministic.m',...
                'ug_solvers_eq.m',...
                'parallel.m',...
                'ug_options.m',...
                'recs_functions.m',...
                'demos.m',...
                'pathnotinstalled.m'};
addpath(fullfile(recsdirectory,'demos'))
for helpfile=HelpFileList
  if ~strcmp(helpfile{1},'simulate.m')
    publish(helpfile{1},PublishOptions);
  else
    publish(helpfile{1},PublishOptionsNoShow);
  end
end
rmpath(fullfile(recsdirectory,'demos'))

%% Functions
if website
  copyfile(fullfile(recsdirectory,'html','layoutfunctionhelp.css'),...
           fullfile(targetdirectory,'layoutfunctionhelp.css'));
  copyfile(fullfile(recsdirectory,'html','bullet-recs.gif'),...
           fullfile(targetdirectory,'bullet-recs.gif'));
  htmlfilelist = ls(fullfile(targetdirectory,'*.html'));
  for i=1:size(htmlfilelist,1)
    txt = fileread(fullfile(targetdirectory,htmlfilelist(i,:)));
    % Replace matlab:doc
    pattern = '<a href="matlab:isfilepresent\(''fsolve.*tt.fsolve..tt...a>';
    txt = regexprep(txt,pattern,'<tt>fsolve</tt>');
    pattern = '<a href="matlab:doc\(''ncpsolve''\)"..tt.ncpsolve..tt...a>';
    txt = regexprep(txt,pattern,'<tt>ncpsolve</tt>');
    pattern = '<a href="matlab:doc\(''ncpsolve''\)"..tt.ncpsolve..tt...a>';
    txt = regexprep(txt,pattern,'<tt>ncpsolve</tt>');
    pattern = '"matlab:doc\(''(\w*)''\)"';
    txt = regexprep(txt,pattern,'"$1.html"');
    fid = fopen(fullfile(targetdirectory,htmlfilelist(i,:)),'w');
    fprintf(fid,'%s',txt);
    fclose(fid);
  end

  FunctionList = {'recsAccuracy',...
                  'recsAuxiliary',...
                  'recsCheck',...
                  'recsConvert',...
                  'recsDecisionRules',...
                  'recsdemos',...
                  'recsFirstGuess',...
                  'recsinterpinit',...
                  'recsmodel',...
                  'recsSimul',...
                  'recsSolveDeterministicPb',...
                  'recsSolveLocal',...
                  'recsSolveREE',...
                  'recsSolveREEFiniteHorizon',...
                  'recsSS',...
                  'lmmcp',...
                  'nsoli',...,
                  'recspathmcp',...
                  'SA',...
                  'SCP'};
  for fn = FunctionList
    if ~strcmp(fn{1},'recsmodel'),  FunctionName = fn{1};
    else,                           FunctionName = 'recsmodel.recsmodel';
    end
    copyfile(fullfile(recsdirectory,[fn{1} '.m']),fullfile(targetdirectory,[fn{1} '.m']));
    txt = help2html(FunctionName);
    % Replace see also
    pattern = '"matlab:helpwin (\w*)"';
    txt = regexprep(txt,pattern,'"$1.html"');
    % Replace recsmodel/recsmodel
    pattern = 'recsmodel/recsmodel';
    txt = strrep(txt,pattern,'recsmodel');
    % Replace view code
    pattern = '"matlab:edit (\w*)"';
    txt = regexprep(txt,pattern,[fn{1} '.m']);
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
stream = RandStream('mt19937ar');
RandStream.setGlobalStream(stream);

currentfolder = cd(fullfile(recsdirectory,'demos'));
DemoFileList = {'cs1','cs2','gro1','gro2','gro3',...
                'sto1','sto2','sto3','sto4','sto5','sto6'};
for demo=DemoFileList
  publish('clearpublish.m',PublishOptions);
  reset(stream);
  publish([demo{1} '.m'],PublishOptions);
end
YamlFileList = {'cs1','gro1','gro2','gro3','sto1','sto2','sto4','sto5','sto6'};
for yaml=YamlFileList
  copyfile([yaml{1} '.yaml'],fullfile(targetdirectory,[yaml{1} '.txt']));
  publish([yaml{1} 'model.m'],PublishOptionsNoExec);
end
delete(fullfile(targetdirectory,'clearpublish.html'));
cd(currentfolder)

%% Build search database
if ~website
  if exist(fullfile(targetdirectory,'helpsearch'),'dir')
    rmdir(fullfile(targetdirectory,'helpsearch'),'s')
  end
  builddocsearchdb(targetdirectory);
end

close all