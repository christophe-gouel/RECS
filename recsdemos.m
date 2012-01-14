% RECSDEMOS runs all RECS demonstration files

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

warning('off','backtrace')

recsdirectory = fileparts(which('recsSimul'));

addpath(fullfile(recsdirectory,'demos'))

disp('GRO1 Stochastic growth model')
clear interp model options
gro1

disp('CS1 Consumption/Saving model with borrowing constraint');
clear interp model options
cs1

disp('STO1 Competitive storage model');
clear interp model options
sto1

disp('STO2 Competitive storage with floor-price backed by public storage');
clear interp model options
sto2

disp('STO3 Competitive storage with price-band backed by public storage');
clear interp model options
sto3

disp('STO4 One small-country storage-trade model');
clear interp model options
sto4

disp('STO5 Two-country storage-trade model');
clear interp model options
sto5

rmpath(fullfile(recsdirectory,'demos'))
