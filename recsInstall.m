function recsInstall
% RECSINSTALL Finalizes RECS installation
%
% RECSINSTALL does three things:
%   - it copies mex files from CompEcon, so the installation of
%     CompEcon must be complete before running RECSINSTALL;
%   - it prepares the html help files;
%   - it download a dolo-recs executable from google code
%     (http://code.google.com/p/dynare-python/downloads/list). Only
%     the executable corresponding to the platform in use is
%     downloaded, so if you use several OS, you need to launch
%     RECSINSTALL on all of them.
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialize
recsdirectory     = fileparts(which('recsSimul'));
compecondirectory = fileparts(which('remsolve'));

%% Copy CompEcon mex files
s1 = copyfile(fullfile(compecondirectory,'private','arraymult.mex*'),recsdirectory);
s2 = copyfile(fullfile(compecondirectory,'arraymult.mex*'),recsdirectory);

if ~s1 && ~s2
  warning(['Failure to copy CompEcon mex files. RECS is not properly ' ...
           'installed, see installation instructions.']);
end

%% Publish html help files
addpath(fullfile(recsdirectory,'html'))
Publish_recs_help

%% Install binary files
if ispc
  extension = '.exe';
else
  extension = '';
end
[~,~] = mkdir(recsdirectory,fullfile('exe',computer('arch')));
urlwrite(['http://dynare-python.googlecode.com/files/dolo-recs-' computer('arch') extension],...
         fullfile(recsdirectory,'exe',computer('arch'),['dolo-recs' ...
                    extension]));