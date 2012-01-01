function filetohelp(filename)
% FILETOHELP

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory = fileparts(which('recsSimul'));
addpath(fullfile(recsdirectory,'html'))

web(['file://' which(filename)],'-helpbrowser');
