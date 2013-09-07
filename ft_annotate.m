function dataout = ft_annotate(cfg, datain)

% FT_ANNOTATE returns the same output data as the user has provided as input, but allows
% to add comments to that data structure. These comments are stored along with the other
% provenance information and can be displayed with FT_ANALYSISPIPELINE. Adding comments
% is especially useful if you have manually (i.e. in plain MATLAB) modified ythe data
% structure, whereby some provenance information is missing.
%
% Use as
%   outdata = ft_examplefunction(cfg, indata)
% where the input data structure can be any of the FieldTrip data structures and where
% cfg is a configuratioun structure that should contain
%
%  cfg.comment    = string
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_ANALYSISPIPELINE, FT_MATH

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


revision = '$Id$';

% do the general setup of the function
ft_defaults                 % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init            % this will show the function help if nargin==0 and return an error
ft_preamble provenance      % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig     % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug           % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar datain  % this reads the input data in case the user specified the cfg.inputfile option

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', 'comment');

% simply copy the input data to the output data
dataout = datain;


% this line is meant to provide some information
% but also to ensure that trackconfig does not remove the cfg.comment field
fprintf('adding the comment: %s\n', cfg.comment);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug            % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig      % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance       % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous datain  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
