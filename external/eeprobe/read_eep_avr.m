function [dat] = read_eep_avr(fn);

% READ_EEP_AVR reads averaged EEG data from an EEProbe *.avr file
% and returns a structure containing the header and data information.
%
% eeg = read_eep_avr(filename)
%
% eeg.label     ... array of labels (1 x nchan)
% eeg.rate      ... sample rate (Hz)
% eeg.npnt      ... number of data points
% eeg.nchan     ... number of channels
% eeg.nsweeps   ... number of trials averaged
% eeg.xmin      ... 
% eeg.xmax      ... 
% eeg.time      ... array of time (1 x npnt)
% eeg.data      ... data array (nchan x npnt)
% eeg.variance  ... variance (nchan x npnt)
% eeg.condlab   ... string with condition label
% eeg.condcol   ... string with color code for condition
% eeg.psi       ... pre-stimulus interval
% eeg.trialc    ... total number of trial in original data
% eeg.rejtrialc ... number of rejected trials
%
% Use plot(eeg.time,eeg.data) to plot the traces at all channels
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_CNT, READ_EEP_TRG, READ_EEP_REJ
%

% Copyright (C) 2002, Robert Oostenveld
%                     Aalborg University, Denmark
%                     http://www.smi.auc.dk/~roberto/
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.1  2004/11/19 14:55:38  jwiskerke
% Added files for use with matlab
%
% Revision 1.2  2003/10/24 13:34:41  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% Revision 1.1.1.1  2003/03/11 15:24:51  roberto
% updated help and copyrights
%
% ANT Software BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

error('could not locate mex file');
