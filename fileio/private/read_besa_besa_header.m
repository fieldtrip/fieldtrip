
function [header] = read_besa_besa_header(fname)
%% Reads BESA .besa format header information and skips data
% See formatting document <a href="matlab:web(http://www.besa.de/downloads/file-formats/)">here</a>
% 
% [alldata,file_info,channel_info,tags,events] = readbesa(fname)
% 
% inputs:
%  fname [string] - path to .besa file
% 
% outputs:
%  header [structure] - Header information
% 
% 
% 
% 2015 - Kristopher Anderson, Knight Lab, Helen Wills Neuroscience Institute, University of California, Berkeley

% For debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
warning on;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Open file
[fid,msg] = fopen(fname,'r');
assert(fid~=-1,'ReadBesaMatlab:ErrorOpeningFile',msg);

% Get length of file
fseek(fid,0,'eof');
file_length = ftell(fid);
fseek(fid,0,'bof');

%% Header Block
[~,ofst_BCF1] = read_tag_offset_pair(fid,'BCF1');

% Read data in header block
while ~feof(fid) && ftell(fid) < (8+ofst_BCF1) % 8 for header tag ('BCF1') and header offset (uint32)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'VERS'
      % File version
      header.orig.file_info.besa_file_version = read_chars(fid,current_length);
    case 'OFFM'
      % Index of first 'file main info' block (BFMI)
      BFMI_offset = fread(fid,1,'*int64');
    case 'OFTL'
      % Index of first 'tag list' block (BTAG)
      BTAG_offset = fread(fid,1,'*int64');
    case 'OFBI'
      % Index of first 'channel and location' block (BCAL)
      BCAL_offset = fread(fid,1,'*int64');
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in header block [BCF1]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
  
end

% Check for necessary header data
if ~exist('BFMI_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBFMI','No BFMI block found in header');
end
if ~exist('BTAG_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBTAG','No BTAG block found in header');
end
if ~exist('BCAL_offset','var')
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoHeaderBCAL','No BCAL block found in header');
end

%% 'tag list' blocks
header.orig.tags.next_BTAG_ofst = BTAG_offset;
header.orig.tags.offsets = [];
header.orig.tags.n_tags = 0;
% Keep reading until no more BTAG blocks
while header.orig.tags.next_BTAG_ofst > 0
  header.orig.tags = read_BTAG(fid, file_length, header.orig.tags);
end
header.orig.tags = rmfield(header.orig.tags,'next_BTAG_ofst');

% Check that file is not much shorter than expected
%  This does not take into account length of final block but might still be useful
if(file_length <= header.orig.tags.tags.position(end))
  fclose(fid);
  error('ReadBesaMatlab:ErrorFileTooShort','Expected file at least %d bytes long but file is %d bytes long',header.orig.tags.tags(end).position,file_length);
end

%% 'file main info' blocks
header.orig.file_info.next_BFMI_ofst = BFMI_offset;
header.orig.file_info.offsets = [];
% Keep reading until no more BFMI blocks
while header.orig.file_info.next_BFMI_ofst > 0
  header.orig.file_info = read_BFMI(fid, file_length, header.orig.file_info);
end
header.orig.file_info = rmfield(header.orig.file_info,'next_BFMI_ofst');
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% 'channel and location' blocks
header.orig.channel_info.next_BCAL_ofst = BCAL_offset;
header.orig.channel_info.offsets = [];
% Keep reading until no more BCAL blocks
while header.orig.channel_info.next_BCAL_ofst > 0
  header.orig.channel_info = read_BCAL(fid, file_length, header.orig.channel_info);
end
header.orig.channel_info = rmfield(header.orig.channel_info,'next_BCAL_ofst');
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

if ~isfield(header.orig.channel_info,'n_channels')
  error('ReadBesaMatlab:ErrorNoHeaderNChannels','Missing number of channels in header [BCAL:CHNR]');
end

% Combine info from channel_info.coord_data and channel_info.channel_states to get actual coordinate data
if(isfield(header.orig.channel_info,'channel_states') && isfield(header.orig.channel_info,'coord_data'))
  for channel_n = 1:header.orig.channel_info.n_channels
    %header.orig.channel_info.channel_locations(channel_n) = [];
    header.orig.channel_info.channel_locations(channel_n).x = NaN;
    header.orig.channel_info.channel_locations(channel_n).y = NaN;
    header.orig.channel_info.channel_locations(channel_n).z = NaN;
    header.orig.channel_info.channel_locations(channel_n).xori = NaN; % Orientation
    header.orig.channel_info.channel_locations(channel_n).yori = NaN;
    header.orig.channel_info.channel_locations(channel_n).zori = NaN;
    header.orig.channel_info.channel_locations(channel_n).x2 = NaN; % Second coil
    header.orig.channel_info.channel_locations(channel_n).y2 = NaN;
    header.orig.channel_info.channel_locations(channel_n).z2 = NaN;
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_SCALPELECTRODE || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).x = double(header.orig.channel_info.coord_data(channel_n,1));
      header.orig.channel_info.channel_locations(channel_n).y = double(header.orig.channel_info.coord_data(channel_n,2));
      header.orig.channel_info.channel_locations(channel_n).z = double(header.orig.channel_info.coord_data(channel_n,3));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).xori = double(header.orig.channel_info.coord_data(channel_n,7));
      header.orig.channel_info.channel_locations(channel_n).yori = double(header.orig.channel_info.coord_data(channel_n,8));
      header.orig.channel_info.channel_locations(channel_n).zori = double(header.orig.channel_info.coord_data(channel_n,9));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER || ...
        header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      header.orig.channel_info.channel_locations(channel_n).x2 = double(header.orig.channel_info.coord_data(channel_n,4));
      header.orig.channel_info.channel_locations(channel_n).y2 = double(header.orig.channel_info.coord_data(channel_n,5));
      header.orig.channel_info.channel_locations(channel_n).z2 = double(header.orig.channel_info.coord_data(channel_n,6));
    end
    if( header.orig.channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE )
      if( header.orig.channel_info.channel_locations(channel_n).x2==0 && ...
          header.orig.channel_info.channel_locations(channel_n).y2==0 && ...
          header.orig.channel_info.channel_locations(channel_n).z2==0 )
        header.orig.channel_info.channel_locations(channel_n).x2 = NaN;
        header.orig.channel_info.channel_locations(channel_n).y2 = NaN;
        header.orig.channel_info.channel_locations(channel_n).z2 = NaN;
      end
    end
  end
end

%% Events
% Collect event block info
header.orig.events.offsets = header.orig.tags.tags.position(strcmp(header.orig.tags.tags.type,'BEVT'));
header.orig.events.offsets = sort(header.orig.events.offsets, 'ascend'); % Later blocks overwrite matching events
for block_n = 1:numel(header.orig.events.offsets)
  header.orig.events = read_BEVT(fid, file_length, header.orig.events, header.orig.events.offsets(block_n));
end
% NEED TO IMPLEMENT OVERWRITES %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Reorganize header structure
header.nChans = header.orig.channel_info.n_channels;
if isfield(header.orig.file_info,'s_rate')
  header.Fs = header.orig.file_info.s_rate;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing sample rate in header');
  header.Fs = [];
end
if isfield(header.orig.file_info,'n_samples')
  header.nSamples = header.orig.file_info.n_samples;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing number of samples in header');
  header.nSamples = [];
end

header.nSamplesPre = []; % Continuous data
header.nTrials     = 1;  % Continuous data

%  Channel labels
if isfield(header.orig.channel_info,'channel_labels')
  header.label = header.orig.channel_info.channel_labels;
else
  warning('ReadBesaMatlab:WarningMissingHeaderInfo','Missing channel labels in header.orig. Creating default channel names');
  for channel_n = 1:header.nChans
    header.label{channel_n} = sprintf('chan%03d', channel_n);
  end
end

%  Channel coordinates
if isfield(header.orig.channel_info,'channel_locations')
  for channel_n = 1:header.nChans
    header.elec.label{channel_n} = header.label{channel_n};
    header.elec.pnt(channel_n,1) = header.orig.channel_info.channel_locations(channel_n).x;
    header.elec.pnt(channel_n,2) = header.orig.channel_info.channel_locations(channel_n).y;
    header.elec.pnt(channel_n,3) = header.orig.channel_info.channel_locations(channel_n).z;
  end
end

function tags = read_BTAG(fid, file_length, tags)
%% Read tag block
% tags [structure] - Existing or blank BTAG structure
%                            Blank needs fields:
%                               next_BTAG_ofst - file offset for BTAG to be read
%                               offsets = []
%                               n_tags = 0
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BTAG section
if(fseek(fid,tags.next_BTAG_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BTAG]',tags.next_BTAG_ofst);
end
tags.offsets(end+1) = tags.next_BTAG_ofst;

% Read BTAG tag and offset
[~,tag_block_length] = read_tag_offset_pair(fid,'BTAG');

% Untagged offset to next BTAG section
tags.next_BTAG_ofst = fread(fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(tags.offsets(end))+uint64(tag_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TAGE'
      % Tag list entry
      tags.n_tags = tags.n_tags+1;
      tags.tags.type{tags.n_tags} = fread(fid,4,'*char')';
      tags.tags.position(tags.n_tags) = fread(fid,1,'*uint64');
      tags.tags.n_samples(tags.n_tags) = double(fread(fid,1,'*uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BTAG]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(tag_block_length) + 8; % 8 for tag and offset
if((tags.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from tag block. Should have read %d bytes', ...
    (ftell(fid)-tags.offsets(end))-expected_length,ftell(fid)-tags.offsets(end),expected_length);
end

function file_info = read_BFMI(fid, file_length, file_info)
%% Read file main info block
% file_info [structure] - Existing or blank BFMI structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BFMI to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BFMI section
if(fseek(fid,file_info.next_BFMI_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BFMI]',file_info.next_BFMI_ofst);
end
file_info.offsets(end+1) = file_info.next_BFMI_ofst;

% Read BFMI tag and offset
[~,fileinfo_block_length] = read_tag_offset_pair(fid,'BFMI');

% Untagged offset to next BFMI section
file_info.next_BFMI_ofst = fread(fid,1,'*int64');

% Create staff field if it doesn't exist already. This is necessary because
%   there is no indication of how many staff to expect, so to increment an
%   array, you need an existing array
if(~isfield(file_info,'staff'))
  file_info.staff = [];
end

% Loop through all tags in data section
while ftell(fid) < (uint64(file_info.offsets(end))+uint64(fileinfo_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMT'
      % Total number of samples
      file_info.n_samples = double(fread(fid,1,'*int64'));
    case 'SAMP'
      % Number of samples per second
      file_info.s_rate = fread(fid,1,'*double');
    case 'FINN'
      % Name of the institution
      file_info.institution.name = read_chars(fid,current_length);
    case 'FINA'
      % Address of the institution
      fina_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < fina_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            file_info.institution.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            file_info.institution.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            file_info.institution.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            file_info.institution.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            file_info.institution.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            file_info.institution.phone_number = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:FINA]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear fina_end;
    case 'FENA'
      % Encryption algorithm
      file_info.encryption = read_chars(fid,current_length);
    case 'FCOA'
      % Compression algorithm
      file_info.compression = read_chars(fid,current_length);
    case 'RECD'
      % Recording start date and time
      file_info.recording_date.start = read_chars(fid,current_length);
    case 'RECE'
      % Recording end date and time
      file_info.recording_date.end = read_chars(fid,current_length);
    case 'RECO'
      % Recording offset to GMT
      file_info.recording_date.gmt_offset = fread(fid,1,'*single');
    case 'RECS'
      % Recording system
      file_info.recording_system.name = read_chars(fid,current_length);
    case 'RIBN'
      % Name of the input box
      file_info.recording_system.info = read_chars(fid,current_length);
    case 'RESW'
      % Name of recording software
      file_info.recording_system.software = read_chars(fid,current_length);
    case 'RATC'
      % Amplifier time constant
      file_info.recording_system.time_constant = fread(fid,1,'*single');
    case 'RSEQ'
      % Sequence number
      file_info.sequence_n = double(fread(fid,1,'*uint32'));
    case 'RSID'
      % Session unique identifier
      file_info.session_id = read_chars(fid,current_length);
    case 'RSNR'
      % Session number
      file_info.sequence_n = double(fread(fid,1,'*int32'));
    case 'RSTC'
      % Study comment
      file_info.comment = read_chars(fid,current_length);
    case 'RSTA'
      % Responsible staff
      % This assumes that, for each staff member, all fields are contiguous
      %   Otherwise, the indices may not line up
      file_info.staff(end+1).name = '';
      file_info.staff(end+1).initials = '';
      file_info.staff(end+1).function = '';
      rsta_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < rsta_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'SNAM'
            % Name
            file_info.staff(end).name = read_chars(fid,current_length);
          case 'ASTA'
            % Initials
            file_info.staff(end).initials = read_chars(fid,current_length);
          case 'ACIT'
            % Function
            file_info.staff(end).function = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:RSTA]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear rsta_end;
    case 'PNAF'
      % Subject first name
      file_info.subject.name.first = read_chars(fid,current_length);
    case 'PNAM'
      % Subject middle name
      file_info.subject.name.middle = read_chars(fid,current_length);
    case 'PATN'
      % Subject last name
      file_info.subject.name.last = read_chars(fid,current_length);
    case 'PNAA'
      % Anonymized subject name
      file_info.subject.anon_name = read_chars(fid,current_length);
    case 'PNAT'
      % Subject title
      file_info.subject.title = read_chars(fid,current_length);
    case 'PATD'
      % Subject date of birth
      file_info.subject.birthdate = read_chars(fid,current_length);
    case 'PDOD'
      % Subject date of death
      file_info.subject.deathdate = read_chars(fid,current_length);
    case 'PAGE'
      % Subject gender
      file_info.subject.gender = read_chars(fid,current_length);
    case 'PAWE'
      % Subject weight
      file_info.subject.weight = fread(fid,1,'*single');
    case 'PAHE'
      % Subject height
      file_info.subject.height = fread(fid,1,'*single');
    case 'PAMS'
      % Subject marital status
      file_info.subject.marital_status = read_chars(fid,current_length);
    case 'PAAD'
      % Subject address
      paad_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < paad_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            file_info.subject.address.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            file_info.subject.address.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            file_info.subject.address.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            file_info.subject.address.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            file_info.subject.address.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            file_info.subject.address.phone_number = read_chars(fid,current_length);
          otherwise
            % Unrecognzed tag. Try to skip forward by offset
            warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
            if((ftell(fid)+current_length) <= file_length)
              if(fseek(fid,current_length,'cof') == -1)
                fclose(fid);
                error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI:PAAD]))',current_length);
              end
            else
              fclose(fid);
              error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
            end
        end
      end
      clear paad_end;
    case 'PALA'
      % Subject language
      file_info.subject.language = read_chars(fid,current_length);
    case 'PAMH'
      % Subject medical history
      file_info.subject.medical_history = read_chars(fid,current_length);
    case 'PATC'
      % Subject comment
      file_info.subject.comment = read_chars(fid,current_length);
    case 'PATI'
      % Subject ID
      file_info.subject.id = read_chars(fid,current_length);
    case 'INF1'
      % Additional information 1
      file_info.additional_info.inf1 = read_chars(fid,current_length);
    case 'INF2'
      % Additional information 2
      file_info.additional_info.inf2 = read_chars(fid,current_length);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BFMI]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(fileinfo_block_length) + 8; % 8 for tag and offset
if((file_info.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from file info block. Should have read %d bytes', ...
    (ftell(fid)-file_info.offsets(end))-expected_length,ftell(fid)-file_info.offsets(end),expected_length);
end

function channel_info = read_BCAL(fid, file_length, channel_info)
%% Read channel info block
% channel_info [structure] - Existing or blank BCAL structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BCAL to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BCAL section
if(fseek(fid,channel_info.next_BCAL_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL]',channel_info.next_BCAL_ofst);
end
channel_info.offsets(end+1) = channel_info.next_BCAL_ofst;

% Read BCAL tag and offset
[~,channel_block_length] = read_tag_offset_pair(fid,'BCAL');

% Untagged offset to next BCAL section
channel_info.next_BCAL_ofst = fread(fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(channel_info.offsets(end))+uint64(channel_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'CHFL'
      % Channel flag
      channel_info.channel_flags.flag = fread(fid,1,'*uint32');
      channel_info.channel_flags.BSA_ELECTRODE_COORDINATES_FROM_LABELS = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0001')),'uint32'));
      channel_info.channel_flags.BSA_SUPPRESS_SPHERE_TO_ELLIPSOID_TRANSFORMATION = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0002')),'uint32'));
      channel_info.channel_flags.BSA_ELECTRODE_COORDINATES_ON_SPHERE = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0004')),'uint32'));
      channel_info.channel_flags.BSA_ADAPT_SPHERICAL_EEG_TO_MEG_COORDS = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0008')),'uint32'));
      channel_info.channel_flags.BSA_SOURCE_CHANNELS_DERIVED_FROM_MEG = logical(bitand(channel_info.channel_flags.flag,uint32(hex2dec('0010')),'uint32'));
    case 'CHTS'
      % Channel type and states of a channel with the specified index
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_states(channel_n).flag = fread(fid,1,'*uint32');
      channel_info.channel_states(channel_n).BSA_CHANTYPE_UNDEFINED = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_POLYGRAPHIC = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00010000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_TRIGGER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00020000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_CORTICALGRID = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00040000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_INTRACRANIAL = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00080000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_SCALPELECTRODE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00100000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00200000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00400000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('01000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00800000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_NKC_REFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('02000000')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANTYPE_CHANSTATE_BAD = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000001')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_REFERENCE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000002')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_INTERPOLRECORDED = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00000004')),'uint32'));
      channel_info.channel_states(channel_n).BSA_CHANSTATE_INVISIBLE = logical(bitand(channel_info.channel_states(channel_n).flag,uint32(hex2dec('00001000')),'uint32'));
    case 'CHCO'
      % Channel coordinates in mm
      n_channels = current_length / 4 / 9; % Divide by 4 for *single and by 9 for number of elements per channel
      channel_info.coord_data = zeros(n_channels,9,'single');
      for channel_n = 1:n_channels
        channel_info.coord_data(channel_n,:) = fread(fid,9,'*single');
      end
      % More processing done later to obtain actual coordinates
    case 'CHNR'
      % Total number of channels
      channel_info.n_channels = double(fread(fid,1,'*uint16'));
    case 'CHLA'
      % Channel label of a channel with the specified index
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_labels{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHET'
      % Electrode thickness
      channel_info.electrode_thickness = fread(fid,1,'*single');
    case 'CHSI'
      % Spline interpolation smoothing constant
      channel_info.spline_smoothing_constant = fread(fid,1,'*single');
    case 'CHLS'
      % Least significant bits of data
      channel_info.lsbs = double(fread(fid,current_length/4,'*single'));
      % Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead. %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CHSF'
      % Sampling frequency
      channel_info.s_rates = fread(fid,current_length/4,'*double');
    case 'HCMM'
      % Head center in mm
      channel_info.head_center.x = fread(fid,1,'*single');
      channel_info.head_center.y = fread(fid,1,'*single');
      channel_info.head_center.z = fread(fid,1,'*single');
    case 'HRMM'
      % Head radius in mm
      channel_info.head_radius = fread(fid,1,'*single');
    case 'FIDC'
      % Fiducial coordinates in mm
      channel_info.fiducial.nasion.x = fread(fid,1,'*single');
      channel_info.fiducial.nasion.y = fread(fid,1,'*single');
      channel_info.fiducial.nasion.z = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.x = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.y = fread(fid,1,'*single');
      channel_info.fiducial.left_preauricular.z = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.x = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.y = fread(fid,1,'*single');
      channel_info.fiducial.right_preauricular.z = fread(fid,1,'*single');
    case 'HSPN'
      % Total number of head surface points
      channel_info.n_addn_surf_pnts = double(fread(fid,1,'*int16'));
    case 'HSPC'
      % Head surface point coordinates
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.head_surface_points{channel_n}.x = fread(fid,1,'*single');
      channel_info.head_surface_points{channel_n}.y = fread(fid,1,'*single');
      channel_info.head_surface_points{channel_n}.z = fread(fid,1,'*single');
    case 'HSPD'
      % Head surface point labels
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.head_surface_points{channel_n}.label = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHCU'
      % Channel units
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_units{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHFI'
      % Filter information
      offset_chfi = ftell(fid);
      offset_end_chfi = int64(offset_chfi)+int64(current_length);
      channel_n = 0;
      while ftell(fid) < offset_end_chfi
        channel_n=channel_n+1;
        filter_object_offset = ftell(fid);
        filter_object_size = fread(fid,1,'*uint32');
        channel_info.filter_info(channel_n).state.state = fread(fid,1,'*uint32');
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_ACTIVE      = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000001')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_ACTIVE     = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000002')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_ACTIVE          = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000004')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_BAND_PASS_ACTIVE      = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000008')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_PRE_LOWCUTOFF_ACTIVE  = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('01000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_PRE_SUBTRACT_BASELINE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('02000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_FORWARD    = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_ZERO_PHASE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000010')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_TYPE_BACKWARD   = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000020')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_06DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_12DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000100')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_24DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000200')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_LOWCUTOFF_SLOPE_48DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000300')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_ZERO_PHASE = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000300')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_FORWARD    = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00001000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_TYPE_BACKWARD   = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00002000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_12DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00000000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_06DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00010000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_24DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00020000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_HIGHCUTOFF_SLOPE_48DB = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00030000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_REMOVE_2ND_HARMONIC = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00100000')),'uint32'));
        channel_info.filter_info(channel_n).state.FLT_NOTCH_REMOVE_3RD_HARMONIC = logical(bitand(channel_info.filter_info(channel_n).state.state,uint32(hex2dec('00200000')),'uint32'));
        channel_info.filter_info(channel_n).low_cutoff  = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).high_cutoff = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).notch_freq  = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).notch_width = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).pass_freq   = fread(fid,1,'*single');
        channel_info.filter_info(channel_n).pass_width  = fread(fid,1,'*single');
        if(ftell(fid) ~= filter_object_offset+filter_object_size)
          warning('ReadBesaMatlab:WarningFilterInfoObject','Did not read expected number bytes in filter object [BCAL:CHFI]. Filter information may be incorrect');
        end
      end
      % Check that expected amout of file was read and move to correct position if not
      if(ftell(fid) ~= offset_end_chfi)
        warning('ReadBesaMatlab:WarningFilterInfoBlock','Did not read expected number of bytes in filter info block [BCAL:CHFI]. Filter information may be incorrect. Skipping to next block');
        fseek(fid,offset_end_chfi+1,'bof');
      end
      % Somewhat complicated structure, no test data %%%%%%%%%%%%%%%%%% TODO
    case 'CHNU'
      % Channel numbers
      channel_info.channel_ns = fread(fid,current_length/4,'*int32');
    case 'CHCM'
      % Channel comments
      channel_n = double(fread(fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      channel_info.channel_comments{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'COMC'
      % BESA CTF component. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    case 'COMH'
      % BESA head transformation. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    case 'CHSC'
      % BESA spatial components. Internal use only
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL:%s]',current_length,current_tag);
      end
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag [BCAL]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Check that expected amout of file was read
expected_length = double(channel_block_length) + 8; % 8 for tag and offset
if((channel_info.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from channel block. Should have read %d bytes', ...
    (ftell(fid)-channel_info.offsets(end))-expected_length,ftell(fid)-channel_info.offsets(end),expected_length);
end


%% EVENT BLOCK FUNCTIONS

function events = read_BEVT(fid, file_length, events, BEVT_offset)
%% Read event block
% BEVT_offset [scalar] - offset of current event block start
% events [structure] - Existing or blank BEVT structure
%                            Blank needs fields:
%                               offsets - sorted array of location of all event blocks
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BEVT section
if(fseek(fid,BEVT_offset,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BEVT]',BEVT_offset);
end

% Read BEVT tag and offset
[~,event_block_length] = read_tag_offset_pair(fid,'BEVT');

% Read LIST tag and offset but don't save anything
read_tag_offset_pair(fid,'LIST');

% Now inside of LIST block
% Read HEAD tag - Assuming that it is first tag in LIST block
[~,head_length] = read_tag_offset_pair(fid,'HEAD');
head_offset = ftell(fid);
% Read data in header block
while ~feof(fid) && ftell(fid) < (head_offset+head_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'EVTS'
      events.n_events = double(fread(fid,1,'*uint32'));
    case 'VERS'
      events.version = double(fread(fid,1,'*uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if((ftell(fid)+current_length) <= file_length)
        if(fseek(fid,current_length,'cof') == -1)
          fclose(fid);
          error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:HEAD]))',current_length);
        end
      else
        fclose(fid);
        error('ReadBesaMatlab:ErrorSkippingForwardAfterUnexpectedTag','Offset after unexpected [%d] tag points to beyond eof [%d]',current_length,file_length);
      end
  end
end

% Read all events as structures and put them into a cell array
for event_n = 1:events.n_events
  [current_tag,current_length] = read_tag_offset_pair(fid);
  events.events{event_n} = read_event_tag(fid,current_tag,current_length);
end

% Check that expected amout of file was read
expected_length = double(event_block_length) + 8; % 8 for tag and offset
if((BEVT_offset+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from event block. Should have read %d bytes', ...
    (ftell(fid)-BEVT_offset)-expected_length,ftell(fid)-BEVT_offset,expected_length);
end

function event_obj = read_event_tag(fid,event_tag,event_length)
% Read an event into a structure
% Create the event object
event_obj.evttype = event_tag;

switch event_tag
  case 'BASE'
    % Base event tag
    event_obj = read_event_tag_base(fid,ftell(fid),event_length,event_obj);
  case 'COMM'
    % Comment event tag
    event_obj = read_event_tag_comm(fid,ftell(fid),event_length,event_obj);
  case 'MARK'
    % Marker event tag
    event_obj = read_event_tag_mark(fid,ftell(fid),event_length,event_obj);
  case 'GENE'
    % Generic event tag
    event_obj = read_event_tag_gene(fid,ftell(fid),event_length,event_obj);
  case 'SEGM'
    % Segment event tag
    event_obj = read_event_tag_segm(fid,ftell(fid),event_length,event_obj);
  case 'ASGM'
    % Average segment start event tag
    event_obj = read_event_tag_asgm(fid,ftell(fid),event_length,event_obj);
  case 'MPS'
    % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    % Multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPS]',event_length);
    end
  case 'MPSC'
    % Classified multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPSC]',event_length);
    end
  case 'PATT'
    % Pattern event tag
    event_obj = read_event_tag_patt(fid,ftell(fid),event_length,event_obj);
  case 'TRIG'
    % Trigger event tag
    event_obj = read_event_tag_trig(fid,ftell(fid),event_length,event_obj);
  case 'PAIR'
    % Paired event tag
    event_obj = read_event_tag_pair(fid,ftell(fid),event_length,event_obj);
  case 'ARTI'
    % Artifact event tag
    event_obj = read_event_tag_arti(fid,ftell(fid),event_length,event_obj);
  case 'EPOC'
    % Epoch event tag
    event_obj = read_event_tag_epoc(fid,ftell(fid),event_length,event_obj);
  case 'IMP'
    % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    % Impedance event tag
    event_obj = read_event_tag_imp(fid,ftell(fid),event_length,event_obj);
  otherwise
    % Unrecognzed tag. Try to skip forward by offset
    warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',event_tag,ftell(fid));
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [LIST]))',event_length);
    end
end

function event_obj = read_event_tag_base(fid,base_offset,base_length, event_obj)
% Read data in the COMM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (base_offset+base_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMP'
      % Sample index (zero based)
      event_obj.sample_n = fread(fid,1,'*int64');
    case 'TIME'
      % Event time
      event_obj.time.year = fread(fid,1,'*int16');
      event_obj.time.month = fread(fid,1,'*int16');
      event_obj.time.dayOfWeek = fread(fid,1,'*int16');
      event_obj.time.day = fread(fid,1,'*int16');
      event_obj.time.hour = fread(fid,1,'*int16');
      event_obj.time.minute = fread(fid,1,'*int16');
      event_obj.time.second = fread(fid,1,'*int16');
      event_obj.time.milliseconds = fread(fid,1,'*int16');
      event_obj.time.microseconds = fread(fid,1,'*single');
      event_obj.time.stateFlag = fread(fid,1,'*uint64');
    case 'SIDX'
      % Segment index (zero based)
      event_obj.segment_index = fread(fid,1,'*int32');
    case 'CODE'
      % Event code (zero based)
      % This value is used by events of type Pattern (PATT) and Trigger (TRIG)
      %  to store the pattern number and the trigger code. Note that the number/code minus 1 is stored.
      %  Additionally, events of type Artifact (ARTI) and Epoch (EPOC) use the code value
      %  internally due to historical reasons. Other event types may use the code value to
      %  store additional information.
      event_obj.code = fread(fid,1,'*int32');
    case 'EVID'
      % Internal BESA event ID (zero based)
      event_obj.besa_event_id = fread(fid,1,'*int32');
    case 'STAT'
      % Event state
      event_obj.state.value = fread(fid,1,'*uint32');
      event_obj.state.EVT_STATE_MARKED1 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000010')),'uint32'));
      event_obj.state.EVT_STATE_MARKED2 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000020')),'uint32'));
      event_obj.state.EVT_STATE_MARKED3 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000040')),'uint32'));
      event_obj.state.EVT_STATE_MARKED4 = logical(bitand(event_obj.state.value,uint32(hex2dec('00000080')),'uint32'));
      event_obj.state.EVT_STATE_DELETED = logical(bitand(event_obj.state.value,uint32(hex2dec('01000000')),'uint32'));
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:BASE]))',current_length);
      end
  end
end

function event_obj = read_event_tag_comm(fid,comm_offset,comm_length, event_obj)
% Read data in the COMM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (comm_offset+comm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TEXT'
      % Event text
      event_obj.text = read_chars(fid,current_length);
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:COMM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_mark(fid,mark_offset,mark_length, event_obj)
% Read data in the MARK event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (mark_offset+mark_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:MARK]))',current_length);
      end
  end
end

function event_obj = read_event_tag_gene(fid,gene_offset,gene_length, event_obj)
% Read data in the GENE event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (gene_offset+gene_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:GENE]))',current_length);
      end
  end
end

function event_obj = read_event_tag_segm(fid,segm_offset,segm_length, event_obj)
% Read data in the SEGM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (segm_offset+segm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SBEG'
      % Segment start time
      event_obj.segment_start.year = fread(fid,1,'*int16');
      event_obj.segment_start.month = fread(fid,1,'*int16');
      event_obj.segment_start.dayOfWeek = fread(fid,1,'*int16');
      event_obj.segment_start.day = fread(fid,1,'*int16');
      event_obj.segment_start.hour = fread(fid,1,'*int16');
      event_obj.segment_start.minute = fread(fid,1,'*int16');
      event_obj.segment_start.second = fread(fid,1,'*int16');
      event_obj.segment_start.milliseconds = fread(fid,1,'*int16');
      event_obj.segment_start.microseconds = fread(fid,1,'*single');
      event_obj.segment_start.stateFlag = fread(fid,1,'*uint64');
    case 'DAYT'
      % Day time of segment start in microseconds
      event_obj.segment_start.dayt = fread(fid,1,'*double');
    case 'INTE'
      % Sampling interval in microseconds
      event_obj.segment_start.sampling_interval = fread(fid,1,'*double');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:SEGM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_asgm(fid,asgm_offset,asgm_length, event_obj)
% Read data in the ASGM event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (asgm_offset+asgm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'STIM'
      % Prestimulus baseline interval in microseconds
      event_obj.baseline_interval = fread(fid,1,'*double');
    case 'AVRS'
      % Number of averages
      event_obj.n_averages = fread(fid,1,'*int32');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:ASGM]))',current_length);
      end
  end
end

function event_obj = read_event_tag_patt(fid,patt_offset,patt_length, event_obj)
% Read data in the PATT event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (patt_offset+patt_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PATT]))',current_length);
      end
  end
end

function event_obj = read_event_tag_trig(fid,trig_offset,trig_length, event_obj)
% Read data in the TRIG event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (trig_offset+trig_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'CODE'
      % Event reaction code
      event_obj.reaction_code = fread(fid,1,'*int32');
    case 'TIME'
      % Event reaction time in seconds
      event_obj.reaction_time = fread(fid,1,'*single');
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:TRIG]))',current_length);
      end
  end
end

function event_obj = read_event_tag_pair(fid,pair_offset,pair_length, event_obj)
% Read data in the PAIR event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (pair_offset+pair_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PART'
      % Event information of the partner event.
      %   The data section is used to write the partner event information
      %   as a data element (starting with <eventtype> tag ID of partner event).
      [event_tag,event_length] = read_tag_offset_pair(fid);
      switch event_tag
        case 'BASE'
          event_obj.partner_event = read_event_tag_base(fid,ftell(fid),event_length,event_obj);
        case 'COMM'
          event_obj.partner_event = read_event_tag_comm(fid,ftell(fid),event_length,event_obj);
        case 'MARK'
          event_obj.partner_event = read_event_tag_mark(fid,ftell(fid),event_length,event_obj);
        case 'GENE'
          event_obj.partner_event = read_event_tag_gene(fid,ftell(fid),event_length,event_obj);
        case 'SEGM'
          event_obj.partner_event = read_event_tag_segm(fid,ftell(fid),event_length,event_obj);
        case 'ASGM'
          event_obj.partner_event = read_event_tag_asgm(fid,ftell(fid),event_length,event_obj);
        case 'MPS'
          % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPS]',event_length);
          end
        case 'MPSC'
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPSC]',event_length);
          end
        case 'PATT'
          event_obj.partner_event = read_event_tag_patt(fid,ftell(fid),event_length,event_obj);
        case 'TRIG'
          event_obj.partner_event = read_event_tag_trig(fid,ftell(fid),event_length,event_obj);
        case 'PAIR'
          event_obj.partner_event = read_event_tag_pair(fid,ftell(fid),event_length,event_obj);
        case 'ARTI'
          event_obj.partner_event = read_event_tag_arti(fid,ftell(fid),event_length,event_obj);
        case 'EPOC'
          event_obj.partner_event = read_event_tag_epoc(fid,ftell(fid),event_length,event_obj);
        case 'IMP'
          % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
          event_obj.partner_event = read_event_tag_imp(fid,ftell(fid),event_length,event_obj);
        otherwise
          % Unrecognzed tag. Try to skip forward by offset
          warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] in PAIR:PART at offset %d',event_tag,ftell(fid));
          if(fseek(fid,event_length,'cof') == -1)
            fclose(fid);
            error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PAIR:PART]))',event_length);
          end
      end
    case 'COMM'
      event_obj = read_event_tag_comm(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:PAIR]))',current_length);
      end
  end
end

function event_obj = read_event_tag_arti(fid,arti_offset,arti_length, event_obj)
% Read data in the ARTI event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (arti_offset+arti_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      event_obj = read_event_tag_pair(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:ARTI]))',current_length);
      end
  end
end

function event_obj = read_event_tag_epoc(fid,epoc_offset,epoc_length, event_obj)
% Read data in the EPOC event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (epoc_offset+epoc_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      event_obj = read_event_tag_pair(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:EPOC]))',current_length);
      end
  end
end

function event_obj = read_event_tag_imp(fid,imp_offset,imp_length, event_obj)
% Read data in the IMP event tag block
% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (imp_offset+imp_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'FORM'
      % Indicates the format used to store impedance values in VAL
      %  Set to 0 if the impedance status (valid/invalid) is stored
      %  Set to 1 if impedance values (in kOhm) are stored
      event_obj.impedance.format = fread(fid,1,'*int32');
    case 'NR'
      % TWO CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      % Number of channels for which impedance information is stored.
      %   (Number of elements stored in TYPE, LABL and VAL)
      event_obj.impedance.n_channels = fread(fid,1,'*uint32');
    case 'TYPE'
      % Channel types
      % The flags used for channel type description are the same as used for the CHTS data elements (specified in chapter 2.3).
      %   Note: Channel type and channel label are used for identification of channel
      %   for which impedance information is set. (Compare to data elements CHTS and CHLA, specified in chapter 2.3).
      event_obj.impedance.types = fread(fid,event_obj.impedance.n_channels,'*uint32'); % Assumes that n_channels has been set
    case 'LABL'
      % Channel labels
      %   Note: Channel type and channel label are used for identification of channel for which impedance information is set.
      %   (Compare to data elements CHTS and CHLA as specified in chapter 2.3).
      event_obj.impedance.labels = read_chars(fid,current_length);
    case 'VAL'
      % TWO CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      % Impedance values
      % Depending on format set in FORM either an impedance STATUS (ok/not ok) or an impedance VALUE (in kOhm) is stored
      %   A value of -1 means that the impedance is not set or invalid
      event_obj.impedance.values = fread(fid,event_obj.impedance.n_channels,'*single'); % Assumes that n_channels has been set
    case 'BASE'
      event_obj = read_event_tag_base(fid,ftell(fid),current_length,event_obj);
    otherwise
      % Unrecognzed tag. Try to skip forward by offset
      warning('ReadBesaMatlab:WarningUnexpectedTag','Read unexpected tag [%s] at offset %d',current_tag,ftell(fid));
      if(fseek(fid,current_length,'cof') == -1)
        fclose(fid);
        error('ReadBesaMatlab:ErrorFseek','fseek to %d failed (after unexpected tag in [BEVT:IMP]))',current_length);
      end
  end
end

function [out_tag, out_offset] = read_tag_offset_pair(fid,expected_tag)
% Read 4 bytes and check if they match expected value
out_tag = fread(fid,4,'*char')';
if(nargin>1)
  % Compare tag with expected tag
  if ~strcmp(expected_tag,out_tag)
    curr_offset = ftell(fid);
    fclose(fid);
    error('ReadBesaMatlab:ErrorTagMismatch','Expecting [%s] but read [%s] at offset %d',expected_tag,out_tag,curr_offset);
  end
end
% Read offset value following tag
out_offset = fread(fid,1,'*uint32');

function outchars = read_chars(fid,n_chars)
% Read n_chars characters from file at current position
%   Replace null character (aka char(0) or '\x0') with ''
%   Note transpose after fread
outchars = regexprep(fread(fid,n_chars,'*char')','\x0','');


