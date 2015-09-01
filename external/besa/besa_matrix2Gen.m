function status = besa_matrix2Gen(data, sRate, basename, varargin)
% BESA_MATRIX2GEN saves sensor level data matrix in double precision to   
% BESA generic data format (basename.generic and basename.dat).
% The exported data can be loaded in BESA Research using 'File/Open'. 
% Select the generated *.generic file and specify the apropriate 
% coordinate files (*.ela/*.elp/*.sfp/*.pos).
%
% Parameters:
%     [data]
%         A 2D or 3D data matrix to be saved. Dimension of data are either 
%         [NumberOfChannels x NumberOfSamples] or 
%         [NumberOfChannels x NumberOfTrials x NumberOfSamples]
% 
%     [sRate]
%         Sampling rate
% 
%     [basename]
%         Basename of the output files (*.generic and *.dat). If data has
%         three dimensions, each trial is exported to a separate output 
%         file basename-1.generic, basename-2.generic etc.
% 
%     [preStimTime]
%         Prestimulus time in ms (set to zero if not specified)
%
% Return:
%     [status] 
%         The status of the writing process: 1 if the writing was
%         successful and 0 otherwise.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_matrix2Gen.m
%
% This file is part of MATLAB2BESA.
%
%    MATLAB2BESA is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    MATLAB2BESA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with MATLAB2BESA. If not, see <http://www.gnu.org/licenses/>.
%
% Author: Karsten Hoechstetter
% Created: 2006-02-03

status = 1;

nChan = size(data,1);
fformat = 'double';
        
datadim = ndims(data);

switch datadim
    case 2
        nSamples = size(data,2);
        % Save data
        fid=fopen([basename,'.dat'],'w');
        fwrite(fid,data,fformat);
        fclose(fid);
        % Generate *.generic file
        fid=fopen([basename,'.generic'],'w');
        fprintf(fid,'BESA Generic Data\n');
        fprintf(fid,'nChannels = %i\n',nChan);
        fprintf(fid,'sRate = %f\n',sRate);
        fprintf(fid,'nSamples = %i\n',nSamples);
        fprintf(fid,'format = %s\n',fformat);
        fprintf(fid,'file = %s.dat\n',basename);
        if size(varargin)>0
            fprintf(fid,'Prestimulus = %f\n',varargin{1});
        end
        fclose(fid);
        
    case 3
        nEpochs = size(data,2);
        nSamples = size(data,3)*nEpochs;
        % Save data
        fid=fopen([basename,'.dat'],'w');
        dataneu=squeeze(data(:,1,:));
        for i=2:nEpochs
            dataneu=[dataneu,squeeze(data(:,i,:))];
        end
        fwrite(fid,squeeze(dataneu),fformat);
        fclose(fid);
        % Generate *.generic file
        fid=fopen([basename,'.generic'],'w');
        fprintf(fid,'BESA Generic Data\n');
        fprintf(fid,'nChannels = %i\n',nChan);
        fprintf(fid,'sRate = %f\n',sRate);
        fprintf(fid,'nSamples = %i\n',nSamples);
        fprintf(fid,'format = %s\n',fformat);
        fprintf(fid,'file = %s.dat\n',basename);
        fprintf(fid,'Epochs = %i\n',nEpochs);
        if size(varargin)>0
            fprintf(fid,'Prestimulus = %f\n',varargin{1});
        end
        fclose(fid);

    otherwise
        display('Error: data matrix must have 2 or 3 dimensions');
        status = 0;
        return
end