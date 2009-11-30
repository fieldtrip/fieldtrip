function [EDF]=closedf(EDF)

% EDF=closedf(EDF)
% Opens an EDF File (European Data Format for Biosignals) into MATLAB
% <A HREF="http://www.medfac.leidenuniv.nl/neurology/knf/kemp/edf.htm">About EDF</A> 
%
% EDF   struct of EDF-Header of a EDF-File

%   Version 2.0
%   15.12.1997
%   Copyright (c) 1997 by  Alois Schloegl
%   a.schloegl@ieee.org 
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%    

EDF.FILE.OPEN=0;
fclose(EDF.FILE.FID);
return;
