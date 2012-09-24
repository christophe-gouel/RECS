del *.zip

md recs-archive
md recs-archive\private
md recs-archive\html
md recs-archive\demos
md recs-archive\exe\win32

REM matlab -wait -nosplash -nojvm -noFigureWindows -r html/Publish_recs_help.m

copy /Y . recs-archive
copy /Y private recs-archive\private
copy /Y html recs-archive\html
copy /Y demos recs-archive\demos
copy /Y exe\win32 recs-archive\exe\win32
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

