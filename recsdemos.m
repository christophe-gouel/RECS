% RECSDEMOS runs all RECS demonstration files

% Copyright (C) 2011-2013 Christophe Gouel
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

disp('GRO3 Stochastic growth model with recursive preferences and stochastic volatility')
clear interp model options
gro3

disp('CS1 Consumption/saving model with borrowing constraint');
clear interp model options
cs1

disp('CS2 Finite horizon consumption/saving model with borrowing constraint');
clear interp model options
cs2

disp('STO1 Competitive storage model');
clear interp model options
sto1

disp('STO1 Competitive storage model with explicit formulation');
clear interp model options
sto1FX

disp('STO2 Competitive storage with floor-price backed by public storage');
clear interp model options
sto2

disp('STO3 Anticipated switch to a public storage policy');
clear interp model options
sto3

disp('STO4 Competitive storage with price-band backed by public storage');
clear interp model options
sto4

disp('STO5 Two-country storage-trade model');
clear interp model options
sto5

disp('STO7 Quarterly storage model with annual inelastic supply');
clear interp model options
sto7

disp('STO7SP Quarterly storage model with informational subperiods and annual inelastic supply');
clear interp model options
sto7SP

rmpath(fullfile(recsdirectory,'demos'))

