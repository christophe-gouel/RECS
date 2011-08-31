function recsInstall
% RECSINSTALL Finalizes RECS installation
%
% RECSINSTALL copies mex files from CompEcon, so the installation
% of CompEcon must be complete before running
% RECSINSTALL. RECSINSTALL prepares also the html help files.

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