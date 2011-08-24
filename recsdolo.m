function recsdolo(inputfile,outputfile)
% RECSDOLO Converts a model from a yaml file to a m file
%
% RECSDOLO uses dolo (https://github.com/albop/dynare-python), a python
% preprocessor, to convert the model described in a yaml file to a file readable by
% MATLAB and RECS programs. In the conversion, dolo calculates the analytic
% representation of all partial derivatives.
%
% RECSDOLO(INPUTFILE,OUTPUTFILE) converts a model in an input file indicated by the
% string INPUTFILE from a yaml format to a m file named by the string OUTPUTFILE.
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin<2, outputfile = strrep(inputfile,'.yaml','model.m'); end
  
recsdirectory = strrep(which('recsSimul'),'recsSimul.m','');

if ~strcmp(computer('arch'),'win32')
  error('Not available on this platform')
end

status = system([recsdirectory 'exe/' computer('arch') '/dolo-recs.exe ' ...
                 inputfile ' ' outputfile]);

if status~=0
  error('Failure to create the model file')
end