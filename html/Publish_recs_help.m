% Publish help pages to html
unix('rm *.png');
options = struct('outputDir','.');
publish('recs_product_page.m',options);
publish('getting_started.m',options);
publish('installation.m',options);
publish('def_dsre.m',options);
publish('user_guide.m',options);
publish('ug_solvers_eq.m',options);
publish('recs_functions.m',options);
publish('pathnotinstalled.m',options);
publish('ug_setting_up.m',options);
