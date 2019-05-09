% mff_importinfo - import information from MFF 'info.xml' file
%
% Usage:
%   info = mff_exportsignal(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  info   - Matlab structure containing informations contained in the MFF
%           file.

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function infoMatlab = mff_importinfo(mffFile)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory         = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

infotype = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Info'));
info = mfffactory.openResourceAtURI( fullfile(mffFile, 'info.xml'), infotype);

infoMatlab.version = 0;
infoMatlab.timezone = [];
infoMatlab.recordtimematlab = 0;
if ~isempty(info)
    %try
        if (info.loadResource() == true)
                        
            % The files version number.
            if (info.getMFFVersionPresent() == true)
                fprintf( 'MFF Version: %d\n', info.getMFFVersion());
                infoMatlab.version = info.getMFFVersion(); % Integer
            end
            
            % The recording time.
            timeVal = char(info.getRecordTime());
            fprintf( 'File''s Recording Time: %s\n', timeVal);
            infoMatlab.recordtime = timeVal;
            
            % get the time zone (duplicate code in mff_importevents and mff_importinfo)
            minusSign = find(timeVal == '+');
            if isempty(minusSign)
                minusSign = find(timeVal == '-');
                minusSign = minusSign(end);
            end
            timeZone = timeVal(minusSign(end):end);
            if length(timeZone) > 6
                timeZone =  [];
                disp('Issue with decoding the time zone');
            end
            infoMatlab.timezone = timeZone;
            
            % decode time
            timeVal  = mff_decodetime(timeVal);
            infoMatlab.recordtimematlab = timeVal;
            
        else
            fprintf( 'Error: Could not load Info resource; file might be corrupted.\n');
        end
        
    %catch
    %   error( 'Unknown while decoding info ressource; send us your data file.\n');
    %end
    
else
    error( 'Error: Could not open the Info resource; check path\n');
end

mfffactory.closeResource(info);
