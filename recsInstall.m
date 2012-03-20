function recsInstall
% RECSINSTALL Finalizes RECS installation
%
% RECSINSTALL does three things:
%   - it copies mex files from CompEcon, so the installation of
%     CompEcon must be complete before running RECSINSTALL;
%   - it downloads a dolo-recs executable from google code
%     (http://code.google.com/p/dynare-python/downloads/list). Only
%     the executable corresponding to the platform in use is
%     downloaded, so if you use several OS, you need to launch
%     RECSINSTALL on all of them.
%   - it prepares the html help files;

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialize
warning('off','backtrace')
recsdirectory     = fileparts(which('recsSimul'));
compecondirectory = fileparts(which('remsolve'));
fprintf('RECS installation:\n')

%% Check recs folder location
if ~isempty(strfind(recsdirectory,' '))
  warning('RECS:SpaceinFolderName',...
          ['RECS folder name, or parent folders, includes spaces, ' ...
           'which may create errors. To avoid errors, please ' ...
           'relocate RECS in a folder without space in its name.'])
end

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

%% Install binary files
fprintf(' - Download dolo-recs executable from Internet: ')
if ispc
  extension = '.exe';
else
  extension = '';
end

if ~(ispc || strcmp(computer('arch'),'glnx86'))
  fprintf('Failure.\n');
  warning('RECS:NoExecForThisArch',...
          'Executable not available on this platform.')
else
  if ispc
    txt = 'win32';
  else
    txt = computer('arch');
  end
  [~,~]  = mkdir(recsdirectory,fullfile('exe',txt));
  [~,s3] = urlwrite(['http://dynare-python.googlecode.com/files/dolo-recs-'...
                     txt extension],...
                    fullfile(recsdirectory,'exe',txt,...
                             ['dolo-recs' extension]));
  if ~s3
  fprintf('Failure.\n');
  warning('RECS:FailureDownloadingURL',...
          ['Failure to download dolo-recs executable. RECS is not properly ' ...
           'installed, see <a href="%s">%s</a> to download dolo-recs manually.'],...
          'http://code.google.com/p/dynare-python/downloads/list',...
          'http://code.google.com/p/dynare-python/downloads/list');
  else
    fprintf('Done.\n');
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
