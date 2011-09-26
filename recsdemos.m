% RECSDEMOS runs all RECS demonstration files

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

warning('off','backtrace')

recsdirectory = fileparts(which('recsSimul'));

addpath([recsdirectory '/demos'])

gro1

disp('CS1 Consumption/Saving model with borrowing constraint');
clear interp model options
cs1

sto1
sto2
sto3
sto4
sto5