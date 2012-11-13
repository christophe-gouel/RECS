del *.zip

md recs-archive
md recs-archive\private
md recs-archive\html
md recs-archive\demos
md recs-archive\Python

cd html
matlab -wait -nosplash -r startup;Publish_recs_help;exit
cd ..

Markdown.pl README.md > README.html

copy /Y . recs-archive
copy /Y private recs-archive\private
robocopy html recs-archive\html /E
copy /Y demos recs-archive\demos
robocopy Python recs-archive\Python /E

cd recs-archive
del .gitignore
del .gitattributes
del /S logfile.tmp 
del Archive.bat
del README.md
del recsInstall.m
7z a ..\recs-%1.zip .
cd ..

del /S /Q recs-archive
rd /S /Q recs-archive

