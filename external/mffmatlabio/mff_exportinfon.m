% mff_exportinfo - export MFF 'info.xml' file from EEGLAB structure
%
% Usage:
%   mff_exportinfon(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)

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

function mff_exportinfon(EEG, mffFile)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.MFFFactory;
import com.egi.services.mff.api.MFFResourceType;
import com.egi.services.mff.api.LocalMFFFactoryDelegate;
import com.egi.services.mff.utility.ResourceUnmarshalException;
import com.egi.services.mff.api.Signal;
import com.egi.services.mff.api.SignalBlock;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
warning('on', 'MATLAB:Java:DuplicateClass');

mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory         = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

infoType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_InfoN'));
if mfffactory.createResourceAtURI(fullfile(mffFile, 'info1.xml'), infoType)
    fprintf('Info1.xml file created successfully\n');
else
    fprintf('Info1.xml ressource already exist, overwriting\n');
end
info = mfffactory.openResourceAtURI( fullfile(mffFile, 'info1.xml'), infoType);

if isfield(EEG.etc, 'infon')
    tmp = javaMethod('valueOf', 'com.egi.services.mff.api.InfoN$InfoNFileType', 'kEEG');
    info.setInfoNFileType(tmp);
    
    tmp = javaObject('com.egi.services.mff.api.InfoNFileTypeEEG');
    info.setInfoNFileTypeInformation(tmp);
    if isfield(EEG.etc.infon, 'infoNFileTypeInformation')
        tmpInfo = EEG.etc.infon.infoNFileTypeInformation;
        if ~isempty(tmpInfo.montageName), tmp.setMontageName(tmpInfo.montageName); end
        if ~isempty(tmpInfo.sensorLayoutName), tmp.setSensorLayoutName(tmpInfo.sensorLayoutName); end
        if ~isempty(tmpInfo.referenceScheme), tmp.setReferenceScheme(tmpInfo.referenceScheme); end
    end
end

info.saveResource();
