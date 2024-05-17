function ks25_phy_toBinary( varargin )

% Read drift corrected data from ks25, unwhiten and rescale for input 
% into C_Waves. Requires a full phy folder of results.

% IMPORTANT NOTE: the preprocessed file only contains channels that KS2.5
% uses for sorting. Channels with connected=0 in the chanMap file OR
% eliminated due to low spike count (ops.minfr_goodchannels > 0) will not
% be included. The information in the phy output is used to generate a new 
% SGLX metadata file which includes the correct channels.

% Another important note: The correct information to unwhiten and rescale
% the KS2.5 drift-corrected output is ONLY available in the rez2.mat output
% file. The whitening.mat.npy and whitening_mat_inv.npy are just identity
% matrices. Voltage values will be incorrectly scaled and the spatial
% fields are post-whitening if the rez2.mat file is unavailable.

% Requires: 
% SGLX_readMeta class: https://github.com/jenniferColonell/SpikeGLX_Datafile_Tools/
% C_Waves: https://billkarsh.github.io/SpikeGLX/#post-processing-tools
% npy-matlab: https://github.com/kwikteam/npy-matlab/tree/master/npy-matlab
% 

% user parameters
modStr = 'ksprocUW'; % will be added to the run name of the output binary and metadata
bRunCWaves = 1; % whether to run C_Waves
exePath = 'C:\Users\colonellj\Documents\C_Waves-win'; % path to C_Waves directory on local machine

if isempty(varargin )
    phyDir = uigetdir([],'select phy directory');
    [binName,binPath] = uigetfile({'*.bin';'*.dat'},'Select PROCESSED binary file');
    fprocFullPath = fullfile(binPath,binName);
    [metaName, mPath] = uigetfile('*.meta','Select ORIGINAL metadata file');
    metaName = metaName(1:length(metaName)-5); % remove extension
else
    % called from another function
    inputCell = varargin(1);
    phyDir = inputCell{1};
    inputCell = varargin(2);
    fprocFullPath = inputCell{1};
    inputCell = varargin(3);
    origMetaFullpath = inputCell{1};
    [mPath, metaName, ~] = fileparts(origMetaFullpath);
end


[baseName,gateStr,triggerStr,probeStr,~] = SGLX_readMeta.parseFileName(metaName);

outName = sprintf( '%s_%s_%s_%s.%s', baseName, modStr, gateStr, triggerStr, probeStr);
binOutName = sprintf( '%s%s', outName, '.ap.bin');
metaOutName = sprintf( '%s%s', outName, '.ap.meta');
outFullPath = fullfile(phyDir, binOutName);


% KS2.5 makes a readable binary from its datashifted data. Rather
% than overlapping the batches, it reads in some extra points for filtering
% and then trims them back off. ops.ntbuff = ntb
% read NTbuff = NT + 3*ntb points
% for a standard batch (neither first nor last) read:
%       --ntb points overlapping the last batch
%       --NT points that belong to this batch
%       --2*ntb more points; first set of ntb will be blended with next
%       batch, 2nd ntb just filtering buffer
% After fitering, points ntb+1:ntb are blended with "dat_prev" which is
% NT+ntb+1:NT+2*ntb saved from the previous batch.
% Batch zero gets ntb zeros prepended to its data and blended with the
% initialized dat_prev (also zeros). 
% After filtering, the data is whitened and then transposed back to NChan
% rows by NT columns to save. When these batches are read for sorting in
% learnTemplates, the data is transposed after reading.

% read whitening matrix and chanMap from phy directory
chanMap = readNPY(fullfile(phyDir,'channel_map.npy'));  % 0 based indicies of the channels in the output file
Nchan = numel(chanMap);

if isfile(fullfile(phyDir,'rez2.mat'))
    rez = load(fullfile(phyDir,'rez2.mat'));
    Wrot = rez.rez.Wrot;
else
    fprintf('Wrot unavailable, will use identity matrix\n');
    fprintf('Voltages will NOT be correct.\n')
    Wrot = eye(Nchan);
    % In the release version of KS2.5, the whitening matrix saved in
    % whitening_mat.npy  is just the identity matrix/rez.ops.scaleProc.
    % whitening_mat_inv.npy is the inverse of whitening_mat.
    % Result: whitening_mat.npy is closer to being the inverse of the
    % whitening matrix than whitening_mat_inv.npy, but NEITHER is correct.
    % To obtain correctly scaled, unwhitened waveforms from the drift
    % corrected data, it is necessary to save the rez2.mat file in the
    % script calling KS2.5.
end


% get number of data points
fp = dir(fprocFullPath);
procSizeBytes = fp.bytes;  %double
totNT = procSizeBytes/(Nchan*2); %total number of time points
fprintf('Total time points in binary %d\n', totNT);
if totNT ~= floor(totNT)
    fprintf('binary doesn not match phy output');
    exit;
end

fid = fopen(fprocFullPath, 'r');
fidW = fopen(outFullPath, 'w'); % open for writing processed data, transposed

% NT is set here to the KS default, but the results are independent of the
% value of NT used by KS2.5.
% NOTE: this is not true for KS2.0.
NT = 65000;

Nbatch = floor(totNT/NT);
batchstart = 0:NT:NT*Nbatch; % batches start at these timepoints
for ibatch = 1:Nbatch
    offset = 2 * Nchan*batchstart(ibatch); % binary file offset in bytes
    fseek(fid, offset, 'bof');
    dat = fread(fid, [Nchan NT], '*int16');
    dat = int16(single(dat')/Wrot);
    fwrite(fidW, dat', 'int16'); % write transposed batch to binary, these chunks are in order
end
% get the end of the file
NTlast = totNT - Nbatch*NT;
offset = 2*Nchan*batchstart(Nbatch);
fseek(fid, offset, 'bof');
dat = fread(fid, [Nchan NTlast], '*int16');
dat = int16(single(dat')/Wrot);
fwrite(fidW, dat', 'int16');
fclose(fid);
fclose(fidW);

fp = dir(outFullPath);
newBytes = uint64(fp.bytes);

% get binaary name and path from origMetaFullPath
origBinName = sprintf('%s.bin', metaName);
origMeta = SGLX_readMeta.ReadMeta(origBinName, mPath);

newMeta = origMeta;
newMeta.fileName = sprintf('%s', outFullPath);
newMeta.fileSizeBytes = sprintf('%ld', newBytes);
newMeta.nSavedChans = sprintf('%d', Nchan);   % read from phy channel_map.npy;
newMeta.snsApLfSy = sprintf('%d,0,0', Nchan);

% map channels to original channels. This is required for correct
% intepretation of gains read from the imro table for NP 1.0-like probes

origChans = SGLX_readMeta.OriginalChans(origMeta) - 1; % zero based indicies of channels in the original binary
origChansNew = origChans(chanMap+1); % zero based indicies of channels in the processed binary
newMeta.snsSaveChanSubset = build_sep_str(origChansNew,',');

% snsChanMap, snsGeomMap and snsShankMap all include only the saved
% channels. Need to get the array of entries from the string, index to
% get the channels included in the output file, 
mapTags = {'snsChanMap', 'snsGeomMap', 'snsShankMap'};
for i = 1:numel(mapTags)
    if isfield(origMeta,mapTags{i})
        currMap = origMeta.(mapTags{i});
        mapArr = split(currMap,')');
        mapArr(numel(mapArr)) = []; % removing last entry from final ')'
        % mapArr is original Nchan + 1 long, first entry is the header
        % each entry is '(' plus the string inside the parens. 
        % Entries in the new map are the header (entry 1) + chanMap entries + 2
        new_map_ind = zeros([Nchan+1,1]);
        new_map_ind(1) = 1; % for header entry
        new_map_ind(2:Nchan+1) = chanMap(1:Nchan) + 2;
        mapArrNew = mapArr(new_map_ind);
        if mapTags{i} == 'snsChanMap'
            % chanMap header to match saved entries
            mapArrNew(1) = {sprintf('(%d,0,0)', Nchan)};
        end
        mapNewStr = build_sep_str(mapArrNew,')');
        newMeta.(mapTags{i}) = sprintf('%s)',mapNewStr);  % adding final close paren
    end
end

SGLX_readMeta.writeMeta(newMeta, fullfile(phyDir,metaOutName) );  
fprintf( 'Output file has %d channels\n', Nchan );

if bRunCWaves
    spkFilePath = fullfile(phyDir,'spike_times.npy');
    spkCluPath = fullfile(phyDir,'spike_clusters.npy');
    clusTablePath = fullfile(phyDir,'clus_Table.npy');
    if ~isfile(clusTablePath)
        phy_to_clusTable(phyDir, Wrot);
    end
    call_CWaves( phyDir, outFullPath, spkFilePath, spkCluPath, clusTablePath, exePath );
end

end

function sepStr = build_sep_str(m, sep)
    % for a 1D array m and separator setp, build the string
    ne = numel(m);
    sepStr = '';
    switch class(m)
        case 'cell'
            for i = 1:ne-1
                % get the contents of the cell
                cellElem = m(i);
                cellStr = cellElem{1};
                sepStr = sprintf('%s%s%s',sepStr,cellStr,sep);
            end
            % last element
            cellElem = m(ne);
            cellStr = cellElem{1};
            sepStr = sprintf('%s%s',sepStr,cellStr);
        otherwise
            for i = 1:ne-1
                sepStr = sprintf('%s%g%s',sepStr,m(i),sep);
            end
            % last element
            sepStr = sprintf('%s%g',sepStr,m(ne));
    end
end

function phy_to_clusTable(phyDir, Wrot)
% buld clusTable for C_Waves from information in phy directory
% clus table is an npy file containing:
% 2-col table (uint32): {num_spike, pk-chan}
    templates = readNPY(fullfile(phyDir,'templates.npy'));
    templates_unwh = zeros(size(templates));
    [nUnit,~,~] = size(templates);
    for i = 1:nUnit
        templates_unwh(i,:,:) = squeeze(templates(i,:,:))/Wrot;
    end
    pp_all = squeeze(max(templates_unwh,[],2) - min(templates_unwh,[],2));
    [~,maxChan] = max(pp_all,[],2);
    spikeClusters = readNPY(fullfile(phyDir,'spike_clusters.npy'));
    [counts,labels] = groupcounts(spikeClusters);
    clu_arr = zeros([nUnit,2],'uint32');
    clu_arr(labels+1,1) = uint32(counts);
    clu_arr(:,2) = uint32(maxChan);
    writeNPY(clu_arr, fullfile(phyDir,'clus_Table.npy'))

end

function call_CWaves( inputPath, binPath, spkFilePath, spkCluPath, clusTablePath, exePath )

% build command line to call CWaves

    args = sprintf("-spikeglx_bin=%s", binPath);
    args = sprintf("%s -clus_table_npy=%s", args, clusTablePath);
    args = sprintf("%s -clus_time_npy=%s", args, spkFilePath);
    args = sprintf("%s -clus_lbl_npy=%s", args, spkCluPath);
    args = sprintf("%s -dest=%s", args, inputPath);
    args = sprintf("%s -prefix=ksproc -samples_per_spike=82 -pre_samples=20 -num_spikes=1000 -snr_radius=8 -snr_radius_um=140", args);
    
    cwaves_cmd = sprintf("%s %s", fullfile(exePath,'runit.bat'), args);
    fprintf("%s\n", cwaves_cmd);
    status = system(cwaves_cmd)

% typical command line:    
% -spikeglx_bin=\\dm11\apig\C_waves_test_data\SC024_092319_NP1.0_Midbrain_g0_tcat.imec0.ap.bin ^
% -clus_table_npy=\\dm11\apig\C_waves_test_data\clus_Table.npy ^
% -clus_time_npy=\\dm11\apig\C_waves_test_data\spike_times.npy ^
% -clus_lbl_npy=\\dm11\apig\C_waves_test_data\spike_clusters.npy ^
% -dest=\\dm11\apig\C_waves_test_data\out ^
% -samples_per_spike=82 -pre_samples=20 -num_spikes=1000 -snr_radius=8 -snr_radius_um=140)
end

