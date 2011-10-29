% Publish help pages to html

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
publish('cs1.m',PublishOptions);
publish('sto1.m',PublishOptions);
cd(currentfolder)

builddocsearchdb(fullfile(recsdirectory,'html'));