% RECSDEMOS runs all RECS demonstration files

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

warning('off','backtrace')

recsdirectory = fileparts(which('recsSimul'));

addpath([recsdirectory '/demos'])

disp('GRO1 Stochastic growth model')
clear interp model options
gro1

disp('CS1 Consumption/Saving model with borrowing constraint');
clear interp model options
cs1

disp('STO1 Competitive storage model');
clear interp model options
sto1
sto2
sto3
sto4
sto5