% RECSINSTALL

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

recsdirectory     = strrep(which('recsSimul'),'recsSimul.m','');
compecondirectory = strrep(which('remsolve.m'),'remsolve.m','');

if ispc
  strarray = ['copy ' compecondirectory 'private\arraymult.mex* ' recsdirectory]
  dos(strarray);
  strarray = ['copy ' compecondirectory         'arraymult.mex* ' recsdirectory]
  dos(strarray);
else % unix or mac
  strarray = ['cp ' compecondirectory 'private/arraymult.mex* ' recsdirectory]
  unix(strarray);
  strarray = ['cp ' compecondirectory         'arraymult.mex* ' recsdirectory];
  unix(strarray);
end
