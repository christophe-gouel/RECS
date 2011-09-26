function filetohelp(filename)
% FILETOHELP

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory = fileparts(which('recsSimul'));
addpath(fullfile(recsdirectory,'demos'))

web(['file://' which(filename)],'-helpbrowser');
