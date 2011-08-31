% Publish help pages to html

recsdirectory = fileparts(which('recsSimul'));
targetdirectory = fullfile(recsdirectory,'html');

if ispc, 
  command = 'del ';
else % unix or mac
  command = 'rm ';
end
system([command fullfile(recsdirectory,'html','*.png')]);
options = struct('outputDir',targetdirectory);
publish('recs_product_page.m',options);
publish('getting_started.m',options);
publish('installation.m',options);
publish('def_dsre.m',options);
publish('user_guide.m',options);
publish('ug_solvers_eq.m',options);
publish('recs_functions.m',options);
publish('pathnotinstalled.m',options);
publish('ug_setting_up.m',options);
