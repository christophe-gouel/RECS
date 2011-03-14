% RECSINSTALL

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

dsemdirectory     = strrep(which('recsSimul'),'recsSimul.m','');
compecondirectory = strrep(which('remsolve.m'),'remsolve.m','');

strarray = ['copy ' compecondirectory 'private\arraymult.mex* ' dsemdirectory];
unix(strarray);
strarray = ['copy ' compecondirectory         'arraymult.mex* ' dsemdirectory];
unix(strarray);
