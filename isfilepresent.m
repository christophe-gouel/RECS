function isfilepresent(filename)
% ISFILEPRESENT

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if exist(filename,'file')
  doc(filename);
else
  switch filename
    case 'fsolve'
      web('http://www.mathworks.fr/help/toolbox/optim/ug/fsolve.html','-helpbrowser');
  end
end