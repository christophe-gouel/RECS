% Publish help pages to html

recsdirectory = fileparts(which('recsSimul'));
targetdirectory = fullfile(recsdirectory,'html');

delete(fullfile(recsdirectory,'html','*.png'));
options = struct('outputDir',targetdirectory);
publish('recs_product_page.m',options);
publish('getting_started.m',options);
publish('installation.m',options);
publish('def_sre.m',options);
publish('MCP.m',options);
publish('user_guide.m',options);
publish('ug_solvers_eq.m',options);
publish('recs_functions.m',options);
publish('pathnotinstalled.m',options);
publish('ug_setting_up.m',options);
publish('demos.m',options);

cd('../demos')
publish('cs1.m',options);
cd('../html')

builddocsearchdb(fullfile(recsdirectory,'html'));