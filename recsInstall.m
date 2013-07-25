function recsInstall
% RECSINSTALL Finalizes RECS installation from source
%
% RECSINSTALL does the following things:
%   - Check if recs folder does not include spaces.
%   - Under Windows, it creates a dolo-recs executable.
%   - It prepares the html help files.
%   - Check if there is the solver Path.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialize
warning('off','backtrace')
recsdirectory     = fileparts(which('recsSimul'));
fprintf('RECS installation:\n')

%% Check recs folder location
if ~isempty(strfind(recsdirectory,' '))
  warning('RECS:SpaceinFolderName',...
          ['RECS folder name, or parent folders, includes spaces, ' ...
           'which may create errors. To avoid errors, please ' ...
           'relocate RECS in a folder without space in its name.'])
end

%% Create dolo-recs executable for Windows
fprintf(' - Create dolo-recs executable for Windows: \n')
if ispc
  dolo     = fullfile(recsdirectory,'Python','dolo');
  dolorecs = fullfile(dolo,'bin','dolo-recs');
  status = system(['pyinstaller.py --onefile --out=pyinstaller --paths=' dolo ' ' dolorecs]);
  if ~status
    movefile(fullfile('.','pyinstaller','dist','dolo-recs.exe'),...
             fullfile(dolo,'bin','dolo-recs.exe'),'f');
    rmdir('pyinstaller','s');
    delete('*.log');
  else
    warning('RECS:DoloRecsCreationFailure',...
            'Failure to create dolo-recs executable for Windows')
  end
end

%% Publish html help files
fprintf(' - Create html help: ')
addpath(fullfile(recsdirectory,'html'))
Publish_recs_help
rmpath(fullfile(recsdirectory,'html'))
fprintf('Done.\n');

%% Check if there is the solver Path
disp('Check dependencies:')
if exist('mcppath','file')
  disp(' - Solver PATH is installed.')
  if isequal(getenv('PATH_LICENSE_STRING'),'')
    disp(' - No license, PATH is in demo mode.')
    warning('RECS:PATHDemoMode',...
            ['In demonstration mode, PATH is limited to 300 variables ' ...
             'and 2000 nonzeros, see <a href="%s">%s</a>'],...
            'http://pages.cs.wisc.edu/~ferris/path.html',...
            'http://pages.cs.wisc.edu/~ferris/path.html')
  else
    disp(' - Existing license for PATH.')
  end
else
  disp(' - Solver PATH is not installed.')
end
