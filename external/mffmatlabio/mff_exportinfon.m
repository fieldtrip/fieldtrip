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

function mff_exportinfon(EEG, mffFile, index)

infon = [ 'info' int2str(index) ];

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory         = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

infoType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_InfoN'));
if mfffactory.createResourceAtURI(fullfile(mffFile, [infon  '.xml']), infoType)
    fprintf('%s.xml file created successfully\n', infon);
else
    fprintf('%s.xml ressource already exist, overwriting\n', infon);
end

info = mfffactory.openResourceAtURI( fullfile(mffFile, [infon  '.xml']), infoType);

% set file type 1 or 2
tmpInfo = [];
if isfield(EEG.etc, infon )
    if isfield(EEG.etc.(infon), 'infoNFileTypeInformation')
        tmpInfo = EEG.etc.(infon).infoNFileTypeInformation;
    end
end
if index == 1 && ~isfield(tmpInfo, 'pnsSetName')
    tmp = javaMethod('valueOf', 'com.egi.services.mff.api.InfoN$InfoNFileType', 'kEEG');
    tmp2 = javaObject('com.egi.services.mff.api.InfoNFileTypeEEG');
else
    tmp = javaMethod('valueOf', 'com.egi.services.mff.api.InfoN$InfoNFileType', 'kPNSData');
    tmp2 = javaObject('com.egi.services.mff.api.InfoNFileTypePNSData');
end

info.setInfoNFileType(tmp);
if isfield(EEG.etc, infon )
    if isfield(EEG.etc.(infon), 'infoNFileTypeInformation')
        tmpInfo = EEG.etc.(infon).infoNFileTypeInformation;
        if isfield(tmpInfo, 'montageName'), tmp2.setMontageName(tmpInfo.montageName); end
        if isfield(tmpInfo, 'sensorLayoutName'), tmp2.setSensorLayoutName(tmpInfo.sensorLayoutName); end
        if isfield(tmpInfo, 'referenceScheme'), tmp2.setReferenceScheme(tmpInfo.referenceScheme); end
        if isfield(tmpInfo, 'pnsSetName'), tmp2.setPNSSetName(tmpInfo.pnsSetName); end
    else
        if index == 1
            tmp2.setMontageName('EEGLAB exported montage');
        else
            tmp2.setPNSSetName('EEGLAB exported PNS channels');
        end
    end
else
    tmp2.setMontageName('EEGLAB exported montage');
end
info.setInfoNFileTypeInformation(tmp2);
info.saveResource();
