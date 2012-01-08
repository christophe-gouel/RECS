function filetohelp(filename)
% FILETOHELP displays in Matlab browsers an ascii file

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory = fileparts(which('recsSimul'));
addpath(fullfile(recsdirectory,'html'))

web(['file://' strrep(which(filename),'\','/')],'-helpbrowser');

rmpath(fullfile(recsdirectory,'html'))
