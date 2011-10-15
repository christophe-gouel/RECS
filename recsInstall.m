function recsInstall
% RECSINSTALL Finalizes RECS installation
%
% RECSINSTALL does three things:
%   - it copies mex files from CompEcon, so the installation of
%     CompEcon must be complete before running RECSINSTALL;
%   - it prepares the html help files;
%   - it downloads a dolo-recs executable from google code
%     (http://code.google.com/p/dynare-python/downloads/list). Only
%     the executable corresponding to the platform in use is
%     downloaded, so if you use several OS, you need to launch
%     RECSINSTALL on all of them.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialize
warning('off','backtrace')
recsdirectory     = fileparts(which('recsSimul'));
compecondirectory = fileparts(which('remsolve'));
fprintf('RECS installation:\n')

%% Copy CompEcon mex files
s1 = copyfile(fullfile(compecondirectory,'private','arraymult.mex*'),recsdirectory);
s2 = copyfile(fullfile(compecondirectory,'arraymult.mex*'),recsdirectory);

fprintf(' - Copy arraymult mex file from CompEcon directory: ')
if ~s1 && ~s2
  fprintf('Failure.\n');
  warning('RECS:FailureCopyArraymuly',...
          ['Failure to copy CompEcon mex files. RECS is not properly ' ...
           'installed, see installation instructions.']);
else
  fprintf('Done.\n');
end


%% Publish html help files
fprintf(' - Create html help: ')
addpath(fullfile(recsdirectory,'html'))
Publish_recs_help
fprintf('Done.\n');

%% Install binary files
fprintf(' - Download dolo-recs executable from Internet: ')
if ispc
  extension = '.exe';
else
  extension = '';
end

if ~(strcmp(computer('arch'),'win32') || strcmp(computer('arch'),'glnx86') ||...
     strcmp(computer('arch'),'glnxa64'))
  fprintf('Failure.\n');
  warning('RECS:NoExecForThisArch',...
          'Executable not available on this platform.')
else
  [~,~] = mkdir(recsdirectory,fullfile('exe',computer('arch')));
  urlwrite(['http://dynare-python.googlecode.com/files/dolo-recs-' computer('arch') extension],...
           fullfile(recsdirectory,'exe',computer('arch'),['dolo-recs' ...
                      extension]));
  fprintf('Done.\n');
end

%% Check if there is the solver Path
disp('Check dependencies:')
if exist('pathmcp','file')
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
