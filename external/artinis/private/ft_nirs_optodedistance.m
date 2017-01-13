function optodedistance = ft_nirs_optodedistance(datain)
% FT_NIRS_OPTODEDISTANCE computes distances between pairs of optodes.
%
% Use as
%   distance = ft_nirs_optodedistance(indata)
% where indata is nirs data.
%
% See also FT_NIRS_REFERENCECHANNELSUBTRACTION

% You are using the FieldTrip NIRS toolbox developed and maintained by 
% Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.fieldtriptoolbox.org
% 
% This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% 
% Creative Commons Attribution-ShareAlike 4.0 International License:
% -----------------------------------
% You are free to:
% 
%     Share — copy and redistribute the material in any medium or format
%     Adapt — remix, transform, and build upon the material
%     for any purpose, even commercially.
% 
%     The licensor cannot revoke these freedoms as long as you follow the 
%     license terms.
% 
% Under the following terms:
% 
%     Attribution — You must give appropriate credit, provide a link to 
%                    the license, and indicate if changes were made. You 
%                    may do so in any reasonable manner, but not in any way 
%                    that suggests the licensor endorses you or your use.
% 
%     ShareAlike — If you remix, transform, or build upon the material, 
%                   you must distribute your contributions under the same 
%                   license as the original.
% 
%     No additional restrictions — You may not apply legal terms or 
%                                   technological measures that legally 
%                                   restrict others from doing anything the 
%                                   license permits.
% 
% -----------------------------------
% 
% This toolbox is not to be used for medical or clinical purposes.
% 
% Copyright (c) 2016 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer: 
% Marc van Wanrooij, DCN, http://www.neural-code.com
% Jörn M. Horschig, Artinis Medical Systems BV, http://www.artinis.com
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure that the input data is raw NIRS-data, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', 'raw', 'senstype', 'nirs');

%% Relevant parameters
label	= datain.label; % transformed channel label, combination Receiver and Transmitter
flabel	= datain.opto.fiberlabel; % fiber label
fpos	= datain.opto.fiberpos; % fiber positions

npos	= numel(label);

xf		= fpos(:,1);
yf		= fpos(:,2);

%% determine distance between Receiver and Transmitter fibers
optodedistance				= NaN(npos,1);

for posIdx		= 1:npos
	str			= label{posIdx};  
  c = textscan(str, '%s%s%s', 'Delimiter', {'-', ' '});
  chanRstr = c{1};
  chanTstr = c{2};
  
	idxR		= match_str(flabel,chanRstr); 
	idxT		= match_str(flabel,chanTstr);
	optodedistance(posIdx)	= sqrt( (xf(idxR)-xf(idxT)).^2+(yf(idxR)-yf(idxT)).^2 ); % Pythagorean theorem
end