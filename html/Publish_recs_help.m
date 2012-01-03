% Publish help pages to html

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory   = fileparts(which('recsSimul'));
targetdirectory = fullfile(recsdirectory,'html');

delete(fullfile(recsdirectory,'html','*.png'));
PublishOptions = struct('outputDir',targetdirectory);
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

currentfolder = cd(fullfile(recsdirectory,'demos'));
copyfile('cs1.yaml',fullfile(targetdirectory,'cs1.txt'));
copyfile('cs1model.m',fullfile(targetdirectory,'cs1model.txt'));
copyfile('gro1.yaml',fullfile(targetdirectory,'gro1.txt'));
copyfile('gro1model.m',fullfile(targetdirectory,'gro1model.txt'));
copyfile('sto1model.m',fullfile(targetdirectory,'sto1model.txt'));
publish('cs1.m',PublishOptions);
publish('gro1.m',PublishOptions);
publish('sto1.m',PublishOptions);
cd(currentfolder)

builddocsearchdb(fullfile(recsdirectory,'html'));