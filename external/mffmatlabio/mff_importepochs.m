% mff_importepochs - import information from MFF 'epochs.xml' file
%
% Usage:
%   epochs = mff_importepochs(mffFile, version);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%  vesion  - file version (optional - default is 3)
%
% Output:
%  epochs - Matlab structure containing informations contained in the MFF file.

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

function continuous = mff_importepochs(mffFile, version)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

if nargin < 2
    version = 3;
end
if version == 0
    divider = 1000;
else
    divider = 1;
end

% Create a factory.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');

%% create Segment to load time
epochsRType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Epochs'));
catURI = [ mffFile filesep 'epochs.xml' ];
if ~exist(catURI)
    catURI = [ mffFile filesep 'epoch.xml' ];
end
epochResource = mfffactorydelegate.openResourceAtURI(catURI, epochsRType);
continuous    = [];
if epochResource.loadResource()
    epochs = epochResource.getEpochs();
    fprintf('Importing epoch.xml ressource: %d data segments\n', epochs.size);
    
    for iEpoch = 1:epochs.size
        singleEpoch = epochs.get(iEpoch-1);
        continuous(iEpoch).begintime  = singleEpoch.getBeginTime()/divider;
        continuous(iEpoch).endtime    = singleEpoch.getEndTime()/divider;
        continuous(iEpoch).firstblock = singleEpoch.getFirstBlock();
        continuous(iEpoch).lastblock  = singleEpoch.getLastBlock();
%         fprintf('Continuous portion Begin Time: %d\n', singleEpoch.getBeginTime());
%         fprintf('Continuous portion End Time: %d\n',   singleEpoch.getEndTime());
%         fprintf('Continuous portion Signal Block This Epoch: %d\n', singleEpoch.getFirstBlock());
%         fprintf('Continuous portion Signal Block This Epoch: %d\n', singleEpoch.getLastBlock());
    end
end

