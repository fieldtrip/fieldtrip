% BIOSIG/T200 contains Matlab/Octave functions to access various biosignal dataformats
% For simiplicity we call all supported files "Biosig"-files.  
% For the the list of supported formats see the references below. 
%
%
% A united interface is provided for all data formats:
%  	SOPEN 	opens an Biosignal file (and reads all header information)
%	SREAD	reads data blockwise
%       SEOF	checks end-of-file
%  	STELL	returns position of file handle
%	SSEEK	moves file handle to position
%  	SREWIND moves file handle to beginning 
%	SCLOSE 	closes an biosignal file 
%  	SWRITE 	writes data blocks 
%
%       GETFILETYPE identifies the type (format) of a file. 
%	SLOAD 	Opens, reads and closes signal files. 
%	SSAVE 	Opens, writes and closes signal files. 
%		SLOAD and SSAVE provide a simple interface to signal files. 
%	SAVE2GDF converts data into GDF-format
% 
% UTILITY FUNCTIONS. In general, it is not recommended 
%	to use them directly. Use them only if you absolute sure what 
%	you are doing. You are warned! 
% 	
%       sload('eventcodes.txt')   loads latest version of table for event codes
%	SAVE2BKR
%	OPENLDR
%	BKROPEN
%	CNTOPEN
%	SDFERROR
%	GDFDATATYP
%       PHYSICALUNITS
%       LEADIDCODEXYZ        
%
% REFERENCES: 
% [1] http://pub.ist.ac.at/~schloegl/biosig/
% [2] http://biosig.sf.net/
%

%	$Id$
%	CopyLeft (c) 1997-2006 by Alois Schloegl <a.schloegl@ieee.org>	
%	This is part of the BIOSIG project http://biosig.sf.net/

