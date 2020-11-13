function dataout = ft_annotate(cfg, datain)

% FT_ANNOTATE returns the same output data as the user has provided as input, but allows
% to add comments to that data structure. These comments are stored along with the other
% provenance information and can be displayed with FT_ANALYSISPIPELINE. Adding comments
% is especially useful if you have manually (i.e. in plain MATLAB) modified the data
% structure, whereby some provenance information is missing.
%
% Use as
%   outdata = ft_examplefunction(cfg, indata)
% where the input data structure can be any of the FieldTrip data structures and
% the configuration structure should contain
%   cfg.comment    = string
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
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', 'comment');

% simply copy the input data to the output data
dataout = datain;


% this line is meant to provide some information
% but also to ensure that trackconfig does not remove the cfg.comment field
fprintf('adding the comment: %s\n', cfg.comment);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous datain
ft_postamble provenance dataout
ft_postamble history dataout
ft_postamble savevar dataout
