function [eeg] = read_eep_cnt(fn);

% READ_EEP_CNT reads continuous EEG data from an EEProbe *.cnt file
% and returns a structure containing the header and data information.
%
% eeg = read_eep_cnt(filename, sample1, sample2)
%
% where sample1 and sample2 are the begin and end sample of the data
% to be read.
%
% eeg.label    ... labels of EEG channels
% eeg.rate     ... sampling rate
% eeg.npnt     ... number of sample in data segment
% eeg.nchan    ... number of channels
% eeg.nsample  
% eeg.time     ... array [1 x npnt] of time points (ms)
% eeg.data     ... array [nchan x npnt] containing eeg data (uV) 
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_TRG, READ_EEP_REJ, READ_EEP_AVR
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

% $Log: read_eep_cnt.m,v $
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/19 14:55:38  jwiskerke
% Added files for use with matlab
%
% Revision 1.1  2004/11/19 14:39:29  jwiskerke
% Initial input into cvs. These files should work.
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
% ANT Software BV, The Netherlands, www.ant-software.nl / info@ant-software.nl
%

error('could not locate mex file');
