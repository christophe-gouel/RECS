del *.zip

:: Create dolo-matlab executable for Windows
pyinstaller --onefile --clean --paths=.\Python\dolo .\Python\dolo\bin\dolo-matlab
copy /Y dist\dolo-matlab.exe Python\dolo\bin
rmdir /S /Q dist
rmdir /S /Q build
del *.spec

:: Create empty folders
md recs-archive
md recs-archive\private
md recs-archive\html
md recs-archive\demos
md recs-archive\Python
md recs-archive\recipes

:: Generate RECS documentation
cd html
matlab -wait -nosplash -r addpath('%1');startup;Publish_recs_help;clear;Publish_recs_help(1);exit
cd ..

:: Convert README.md to html
pandoc README.md -s -o README.html

:: Copy all files
copy /Y . recs-archive
copy /Y private recs-archive\private
robocopy html recs-archive\html /E
copy /Y demos recs-archive\demos
robocopy Python\dolo recs-archive\Python\dolo /E
copy /Y recipes recs-archive\recipes

:: Clean the files and compress
cd recs-archive
del /S .git*
del /S logfile.tmp
del /S *.mex*
del Archive.bat
rename README.md README.txt
del recsInstall.m
rmdir /S /Q Python\dolo\.git
7z a ..\recs-%2.zip .
cd ..

:: Clean temporary folder
rmdir /S /Q recs-archive

