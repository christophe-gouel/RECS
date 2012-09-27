function Publish_recs_help
% PUBLISH_RECS_HELP publishes help pages to html

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory   = fileparts(which('recsSimul'));
targetdirectory = fullfile(recsdirectory,'html');

delete(fullfile(recsdirectory,'html','*.png'));
delete(fullfile(recsdirectory,'html','*.txt'));
delete(fullfile(recsdirectory,'html','*.yaml'));
delete(fullfile(recsdirectory,'html','*.html'));
PublishOptions = struct('outputDir',targetdirectory);
PublishOptionsNoExec = struct('outputDir',targetdirectory,...
                              'evalCode',false);

%% Documentation
publish('recs_product_page.m',PublishOptions);
publish('getting_started.m',PublishOptions);
publish('installation.m',PublishOptions);
publish('def_sre.m',PublishOptions);
publish('MCP.m',PublishOptions);
publish('user_guide.m',PublishOptions);
publish('ug_solvers_eq.m',PublishOptions);
publish('recs_functions.m',PublishOptions);
publish('pathnotinstalled.m',PublishOptions);
publish('ug_setting_up.m',PublishOptions);
publish('ug_model_files.m',PublishOptions);
publish('demos.m',PublishOptions);

%% License
copyfile(fullfile(recsdirectory,'LICENSE.txt'),fullfile(targetdirectory,'LICENSE.txt'));

%% Demonstration
currentfolder = cd(fullfile(recsdirectory,'demos'));
copyfile('cs1.yaml',fullfile(targetdirectory,'cs1.yaml'));
copyfile('gro1.yaml',fullfile(targetdirectory,'gro1.yaml'));
copyfile('gro2.yaml',fullfile(targetdirectory,'gro2.yaml'));
publish('cs1model.m',PublishOptionsNoExec);
publish('gro1model.m',PublishOptionsNoExec);
publish('gro2model.m',PublishOptionsNoExec);
publish('sto1model.m',PublishOptionsNoExec);
publish('sto2model.m',PublishOptionsNoExec);
publish('sto4model.m',PublishOptionsNoExec);
publish('sto5model.m',PublishOptionsNoExec);
publish('sto6model.m',PublishOptionsNoExec);
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
publish('sto4.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto5.m',PublishOptions);
publish('clearpublish.m',PublishOptions);
publish('sto6.m',PublishOptions);
delete(fullfile(recsdirectory,'html','clearpublish.html'));
cd(currentfolder)

builddocsearchdb(fullfile(recsdirectory,'html'));