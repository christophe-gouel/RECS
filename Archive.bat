del *.zip

md recs-archive
md recs-archive\private
md recs-archive\html
md recs-archive\demos
md recs-archive\Python

REM matlab -wait -nosplash -nojvm -noFigureWindows -r html/Publish_recs_help.m

copy /Y . recs-archive
copy /Y private recs-archive\private
robocopy html recs-archive\html /E
copy /Y demos recs-archive\demos
robocopy Python recs-archive\Python /E
copy /Y Archive\README.txt recs-archive

cd recs-archive
del .gitignore
del /S logfile.tmp 
del Archive.bat
del README.md
del recsInstall.m
7z a ..\recs-%1.zip .
cd ..

del /S /Q recs-archive
rd /S /Q recs-archive

