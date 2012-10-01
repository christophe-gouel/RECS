% RECSDEMOS runs all RECS demonstration files

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

warning('off','backtrace')

recsdirectory = fileparts(which('recsSimul'));

addpath(fullfile(recsdirectory,'demos'))

disp('GRO1 Stochastic growth model')
clear interp model options
gro1

disp('GRO2 Stochastic growth model with irreversible investment')
clear interp model options
gro2

disp('CS1 Consumption/saving model with borrowing constraint');
clear interp model options
cs1

disp('CS2 Finite horizon consumption/saving model with borrowing constraint');
clear interp model options
cs2

disp('STO1 Competitive storage model');
clear interp model options
sto1

disp('STO2 Competitive storage with floor-price backed by public storage');
clear interp model options
sto2

disp('STO3 Anticipated switch to a public storage policy');
clear interp model options
sto3

disp('STO4 Competitive storage with price-band backed by public storage');
clear interp model options
sto4

disp('STO5 One small-country storage-trade model');
clear interp model options
sto5

disp('STO6 Two-country storage-trade model');
clear interp model options
sto6

rmpath(fullfile(recsdirectory,'demos'))
