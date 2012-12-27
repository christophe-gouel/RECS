del *.zip

:: Create dolo-recs executable for Windows
pyinstaller.py -F --out=pyinstaller --paths=.\Python\dolo .\Python\dolo\bin\dolo-recs.py
copy /Y pyinstaller\dist\dolo-recs.exe Python\dolo\bin
rmdir /S /Q pyinstaller
del *.log

:: Create empty folders
md recs-archive
md recs-archive\private
md recs-archive\html
md recs-archive\demos
md recs-archive\Python

:: Generate RECS documentation
cd html
matlab -wait -nosplash -r startup;Publish_recs_help;exit
cd ..

:: Convert README.md to html
Markdown.pl README.md > README.html

:: Copy all files
copy /Y . recs-archive
copy /Y private recs-archive\private
robocopy html recs-archive\html /E
copy /Y demos recs-archive\demos
robocopy Python recs-archive\Python /E

:: Clean the files and compress
cd recs-archive
del .gitignore
del .gitattributes
del /S logfile.tmp 
del Archive.bat
rename README.md README.txt
del recsInstall.m
del Python\dolo\.gitignore
rmdir /S /Q Python\dolo\.git
7z a ..\recs-%1.zip .
cd ..

:: Clean temporary folder
rmdir /S /Q recs-archive

