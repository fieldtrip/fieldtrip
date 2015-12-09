
function [BFMI_struct,BCAL_struct,BTAG_struct,BEVT_struct,alldata] = readbesa(fname)
%% Reads BESA .besa format files
% See formatting document <a href="matlab:web(http://www.besa.de/downloads/file-formats/)">here</a>
% 
% [alldata,BFMI_struct,BCAL_struct,BTAG_struct,BEVT_struct] = readbesa(fname)
% 
% inputs:
%  fname [string] - path to .besa file
% 
% outputs:
%  alldata [n_samples x n_channels] - The data
%  BFMI_struct [structure] - File info
%  BCAL_struct [structure] - Channel info
%  BTAG_struct [structure] - Tags
%  BEVT_struct [structure] - Events
% TODO: Consolidate these into single header structure
% 
% 
% 
% 2015 - Kristopher Anderson, Knight Lab, Helen Wills Neuroscience Institute, University of California, Berkeley

% For debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
BFMI_struct = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
BCAL_struct = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
BTAG_struct = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
BEVT_struct = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
alldata = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
warning on;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

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
      BFMI_struct.besa_file_version = read_chars(fid,current_length);
    case 'OFFM'
      % Index of first 'file main info' block (BFMI)
      BFMI_offset = freadshow(0, fid,1,'*int64');
    case 'OFTL'
      % Index of first 'tag list' block (BTAG)
      BTAG_offset = freadshow(0, fid,1,'*int64');
    case 'OFBI'
      % Index of first 'channel and location' block (BCAL)
      BCAL_offset = freadshow(0, fid,1,'*int64');
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
BTAG_struct.next_BTAG_ofst = BTAG_offset;
BTAG_struct.offsets = [];
BTAG_struct.n_tags = 0;
% Keep reading until no more BTAG blocks
while BTAG_struct.next_BTAG_ofst > 0
  BTAG_struct = read_BTAG(fid, file_length, BTAG_struct);
end
BTAG_struct = rmfield(BTAG_struct,'next_BTAG_ofst');

% Check that file is not much shorter than expected
%  This does not take into account length of final block but might still be useful
if(file_length <= BTAG_struct.tags.position(end))
  fclose(fid);
  error('ReadBesaMatlab:ErrorFileTooShort','Expected file at least %d bytes long but file is %d bytes long',BTAG_struct.tags(end).position,file_length);
end

%% 'file main info' blocks
BFMI_struct.next_BFMI_ofst = BFMI_offset;
BFMI_struct.offsets = [];
% Keep reading until no more BFMI blocks
while BFMI_struct.next_BFMI_ofst > 0
  BFMI_struct = read_BFMI(fid, file_length, BFMI_struct);
end
BFMI_struct = rmfield(BFMI_struct,'next_BFMI_ofst');
% NEED TO CHECK THAT OVERWRITES ARE HANDLED CORRECTLY %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% 'channel and location' blocks
BCAL_struct.next_BCAL_ofst = BCAL_offset;
BCAL_struct.offsets = [];
% Keep reading until no more BCAL blocks
while BCAL_struct.next_BCAL_ofst > 0
  BCAL_struct = read_BCAL(fid, file_length, BCAL_struct);
end
BCAL_struct = rmfield(BCAL_struct,'next_BCAL_ofst');
% NEED TO CHECK THAT OVERWRITES ARE HANDLED CORRECTLY %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Combine info from BCAL_struct.coord_data and BCAL_struct.channel_states
% to get actual coordinate data
% The location/orientation values are set depending on the channel type that was set in CHTS section for the associated channel:
% 0x00100000 (BSA_CHANTYPE_SCALPELECTRODE) Entries 0 to 2 contain the location (x,y,z).
% Entries 3 to 8 are ignored.
% 0x00200000 (BSA_CHANTYPE_MAGNETOMETER) Entries 0 to 2 contain the location (x,y,z). Entries 3 to 5 are ignored.
% Entries 6 to 8 contain the orientation (x,y,z).
% 0x00400000 (BSA_CHANTYPE_AXIAL_GRADIOMETER) 0x01000000 (BSA_CHANTYPE_PLANAR_GRADIOMETER) Entries 0 to 2 contain the location (x,y,z) of the 1st coil.
% Entries 3 to 5 contain the location of the 2nd coil. Entries 6 to 8 contain the orientation (x,y,z).
% 0x00800000 (BSA_CHANTYPE_MEGREFERENCE) Entries 0 to 2 contain the location (x,y,z). Entries 6 to 8 contain the orientation (x,y,z).
% Entries 3 to 5 are checked to decide if the channel is a magnetometer or gradiometer: If all three coordinates are zero the sensor is assumed to be a magnetometer. Otherwise the sensor is assumed to be a gradiometer and the entries 3 to 5 specify the location of the 2nd coil.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Events
% Collect event block info
BEVT_struct.offsets = BTAG_struct.tags.position(strcmp(BTAG_struct.tags.type,'BEVT'));
BEVT_struct.offsets = sort(BEVT_struct.offsets, 'ascend'); % Later blocks overwrite matching events
for block_n = 1:numel(BEVT_struct.offsets)
  BEVT_struct = read_BEVT(fid, file_length, BEVT_struct, BEVT_struct.offsets(block_n));
  fprintf('Events Block %d of %d\n',block_n,numel(BEVT_struct.offsets)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
end
% NEED TO CHECK THAT OVERWRITES ARE HANDLED CORRECTLY %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

%% Data
% Collect data block info
data_block_offsets   = BTAG_struct.tags.position(strcmp(BTAG_struct.tags.type,'BDAT'));
data_block_n_samples = BTAG_struct.tags.n_samples(strcmp(BTAG_struct.tags.type,'BDAT'));
data_block_samples_begin = ones(numel(data_block_n_samples),1);
data_block_samples_end   = ones(numel(data_block_n_samples),1)*data_block_n_samples(1);
for block_n = 2:numel(data_block_n_samples)
  data_block_samples_begin(block_n) = data_block_samples_end(block_n-1) + 1;
  data_block_samples_end(block_n)   = data_block_samples_end(block_n-1) + data_block_n_samples(block_n);
end

% Check for necessary values
if(~isfield(BCAL_struct,'n_channels'))
  fclose(fid);
  error('ReadBesaMatlab:ErrorNoNChannels','BCAL_struct.n_channels does not exist. This is needed for reading data blocks');
end
if(~isfield(BCAL_struct,'lsbs'))
  % No least significant bit values found, so setting them all to 1.0
  BCAL_struct.lsbs = ones(BCAL_struct.n_channels,1,'double');
end

% Loop over all data blocks and add to alldata matrix
alldata = zeros(sum(data_block_n_samples),BCAL_struct.n_channels);
for block_n = 1:numel(data_block_n_samples)
  alldata(data_block_samples_begin(block_n):data_block_samples_end(block_n),:) = ...
    read_BDAT(fid, file_length, data_block_offsets(block_n), BCAL_struct.n_channels, BCAL_struct.lsbs);
  fprintf('Data Block %d of %d\n',block_n,numel(data_block_n_samples)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
end
% NEED TO CHECK THAT OVERWRITES ARE HANDLED CORRECTLY %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO




function BTAG_struct = read_BTAG(fid, file_length, BTAG_struct)
%% Read tag block
% BTAG_struct [structure] - Existing or blank BTAG structure
%                            Blank needs fields:
%                               next_BTAG_ofst - file offset for BTAG to be read
%                               offsets = []
%                               n_tags = 0
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BTAG section
if(fseek(fid,BTAG_struct.next_BTAG_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BTAG]',BTAG_struct.next_BTAG_ofst);
end
BTAG_struct.offsets(end+1) = BTAG_struct.next_BTAG_ofst;

% Read BTAG tag and offset
[~,tag_block_length] = read_tag_offset_pair(fid,'BTAG');

% Untagged offset to next BTAG section
BTAG_struct.next_BTAG_ofst = freadshow(0, fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(BTAG_struct.offsets(end))+uint64(tag_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TAGE'
      % Tag list entry
      BTAG_struct.n_tags = BTAG_struct.n_tags+1;
      BTAG_struct.tags.type{BTAG_struct.n_tags} = freadshow(0, fid,4,'*char')';
      BTAG_struct.tags.position(BTAG_struct.n_tags) = freadshow(0, fid,1,'*uint64');
      BTAG_struct.tags.n_samples(BTAG_struct.n_tags) = double(freadshow(0, fid,1,'*uint32'));
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
if((BTAG_struct.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from tag block. Should have read %d bytes', ...
    (ftell(fid)-BTAG_struct.offsets(end))-expected_length,ftell(fid)-BTAG_struct.offsets(end),expected_length);
end

function BFMI_struct = read_BFMI(fid, file_length, BFMI_struct)
%% Read file main info block
% BFMI_struct [structure] - Existing or blank BFMI structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BFMI to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BFMI section
if(fseek(fid,BFMI_struct.next_BFMI_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BFMI]',BFMI_struct.next_BFMI_ofst);
end
BFMI_struct.offsets(end+1) = BFMI_struct.next_BFMI_ofst;

% Read BFMI tag and offset
[~,fileinfo_block_length] = read_tag_offset_pair(fid,'BFMI');

% Untagged offset to next BFMI section
BFMI_struct.next_BFMI_ofst = freadshow(0, fid,1,'*int64');

% Create staff field if it doesn't exist already. This is necessary because
%   there is no indication of how many staff to expect, so to increment an
%   array, you need an existing array
if(~isfield(BFMI_struct,'staff'))
  BFMI_struct.staff = [];
end

% Loop through all tags in data section
while ftell(fid) < (uint64(BFMI_struct.offsets(end))+uint64(fileinfo_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMT'
      % Total number of samples
      BFMI_struct.n_samples = double(freadshow(0, fid,1,'*int64'));
    case 'SAMP'
      % Number of samples per second
      BFMI_struct.s_rate = freadshow(0, fid,1,'*double');
    case 'FINN'
      % Name of the institution
      BFMI_struct.institution.name = read_chars(fid,current_length);
    case 'FINA'
      % Address of the institution
      fina_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < fina_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            BFMI_struct.institution.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            BFMI_struct.institution.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            BFMI_struct.institution.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            BFMI_struct.institution.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            BFMI_struct.institution.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            BFMI_struct.institution.phone_number = read_chars(fid,current_length);
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
      BFMI_struct.encryption = read_chars(fid,current_length);
    case 'FCOA'
      % Compression algorithm
      BFMI_struct.compression = read_chars(fid,current_length);
    case 'RECD'
      % Recording start date and time
      BFMI_struct.recording_date.start = read_chars(fid,current_length);
    case 'RECE'
      % Recording end date and time
      BFMI_struct.recording_date.end = read_chars(fid,current_length);
    case 'RECO'
      % Recording offset to GMT
      BFMI_struct.recording_date.gmt_offset = freadshow(0, fid,1,'*single');
    case 'RECS'
      % Recording system
      BFMI_struct.recording_system.name = read_chars(fid,current_length);
    case 'RIBN'
      % Name of the input box
      BFMI_struct.recording_system.info = read_chars(fid,current_length);
    case 'RESW'
      % Name of recording software
      BFMI_struct.recording_system.software = read_chars(fid,current_length);
    case 'RATC'
      % Amplifier time constant
      BFMI_struct.recording_system.time_constant = freadshow(0, fid,1,'*single');
    case 'RSEQ'
      % Sequence number
      BFMI_struct.sequence_n = double(freadshow(0, fid,1,'*uint32'));
    case 'RSID'
      % Session unique identifier
      BFMI_struct.session_id = read_chars(fid,current_length);
    case 'RSNR'
      % Session number
      BFMI_struct.sequence_n = double(freadshow(0, fid,1,'*int32'));
    case 'RSTC'
      % Study comment
      BFMI_struct.comment = read_chars(fid,current_length);
    case 'RSTA'
      % Responsible staff
      % This assumes that, for each staff member, all fields are contiguous
      %   Otherwise, the indices may not line up
      BFMI_struct.staff(end+1).name = '';
      BFMI_struct.staff(end+1).initials = '';
      BFMI_struct.staff(end+1).function = '';
      rsta_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < rsta_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'SNAM'
            % Name
            BFMI_struct.staff(end).name = read_chars(fid,current_length);
          case 'ASTA'
            % Initials
            BFMI_struct.staff(end).initials = read_chars(fid,current_length);
          case 'ACIT'
            % Function
            BFMI_struct.staff(end).function = read_chars(fid,current_length);
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
      BFMI_struct.subject.name.first = read_chars(fid,current_length);
    case 'PNAM'
      % Subject middle name
      BFMI_struct.subject.name.middle = read_chars(fid,current_length);
    case 'PATN'
      % Subject last name
      BFMI_struct.subject.name.last = read_chars(fid,current_length);
    case 'PNAA'
      % Anonymized subject name
      BFMI_struct.subject.anon_name = read_chars(fid,current_length);
    case 'PNAT'
      % Subject title
      BFMI_struct.subject.title = read_chars(fid,current_length);
    case 'PATD'
      % Subject date of birth
      BFMI_struct.subject.birthdate = read_chars(fid,current_length);
    case 'PDOD'
      % Subject date of death
      BFMI_struct.subject.deathdate = read_chars(fid,current_length);
    case 'PAGE'
      % Subject gender
      BFMI_struct.subject.gender = read_chars(fid,current_length);
    case 'PAWE'
      % Subject weight
      BFMI_struct.subject.weight = freadshow(0, fid,1,'*single');
    case 'PAHE'
      % Subject height
      BFMI_struct.subject.height = freadshow(0, fid,1,'*single');
    case 'PAMS'
      % Subject marital status
      BFMI_struct.subject.marital_status = read_chars(fid,current_length);
    case 'PAAD'
      % Subject address
      paad_end = ftell(fid)+current_length;
      while ~feof(fid) && ftell(fid) < paad_end
        [current_tag,current_length] = read_tag_offset_pair(fid);
        switch current_tag
          case 'ASTR'
            % Street name
            BFMI_struct.subject.address.street_name = read_chars(fid,current_length);
          case 'ASTA'
            % State
            BFMI_struct.subject.address.state = read_chars(fid,current_length);
          case 'ACIT'
            % City
            BFMI_struct.subject.address.city = read_chars(fid,current_length);
          case 'APOS'
            % Post code
            BFMI_struct.subject.address.post_code = read_chars(fid,current_length);
          case 'ACOU'
            % Country
            BFMI_struct.subject.address.country = read_chars(fid,current_length);
          case 'APHO'
            % Phone number
            BFMI_struct.subject.address.phone_number = read_chars(fid,current_length);
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
      BFMI_struct.subject.language = read_chars(fid,current_length);
    case 'PAMH'
      % Subject medical history
      BFMI_struct.subject.medical_history = read_chars(fid,current_length);
    case 'PATC'
      % Subject comment
      BFMI_struct.subject.comment = read_chars(fid,current_length);
    case 'PATI'
      % Subject ID
      BFMI_struct.subject.id = read_chars(fid,current_length);
    case 'INF1'
      % Additional information 1
      BFMI_struct.additional_info.inf1 = read_chars(fid,current_length);
    case 'INF2'
      % Additional information 2
      BFMI_struct.additional_info.inf2 = read_chars(fid,current_length);
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
if((BFMI_struct.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from file info block. Should have read %d bytes', ...
    (ftell(fid)-BFMI_struct.offsets(end))-expected_length,ftell(fid)-BFMI_struct.offsets(end),expected_length);
end

function BCAL_struct = read_BCAL(fid, file_length, BCAL_struct)
%% Read channel info block
% BCAL_struct [structure] - Existing or blank BCAL structure
%                            Blank needs fields:
%                               next_BFMI_ofst - file offset for BCAL to be read
%                               offsets = []
% file_length [scalar] - Length of file in bytes
% fid [scalar] - File identifier

% Skip to start of BCAL section
if(fseek(fid,BCAL_struct.next_BCAL_ofst,'bof') == -1)
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BCAL]',BCAL_struct.next_BCAL_ofst);
end
BCAL_struct.offsets(end+1) = BCAL_struct.next_BCAL_ofst;

% Read BCAL tag and offset
[~,channel_block_length] = read_tag_offset_pair(fid,'BCAL');

% Untagged offset to next BCAL section
BCAL_struct.next_BCAL_ofst = freadshow(0, fid,1,'*int64');

% Loop through all tags in data section
while ftell(fid) < (uint64(BCAL_struct.offsets(end))+uint64(channel_block_length))
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'CHFL'
      % Channel flag
      BCAL_struct.channel_flags.flag = freadshow(0, fid,1,'*uint32');
      BCAL_struct.channel_flags.BSA_ELECTRODE_COORDINATES_FROM_LABELS = logical(bitand(BCAL_struct.channel_flags.flag,uint32(hex2dec('0001')),'uint32'));
      BCAL_struct.channel_flags.BSA_SUPPRESS_SPHERE_TO_ELLIPSOID_TRANSFORMATION = logical(bitand(BCAL_struct.channel_flags.flag,uint32(hex2dec('0002')),'uint32'));
      BCAL_struct.channel_flags.BSA_ELECTRODE_COORDINATES_ON_SPHERE = logical(bitand(BCAL_struct.channel_flags.flag,uint32(hex2dec('0004')),'uint32'));
      BCAL_struct.channel_flags.BSA_ADAPT_SPHERICAL_EEG_TO_MEG_COORDS = logical(bitand(BCAL_struct.channel_flags.flag,uint32(hex2dec('0008')),'uint32'));
      BCAL_struct.channel_flags.BSA_SOURCE_CHANNELS_DERIVED_FROM_MEG = logical(bitand(BCAL_struct.channel_flags.flag,uint32(hex2dec('0010')),'uint32'));
    case 'CHTS'
      % Channel type and states of a channel with the specified index
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.channel_states(channel_n).flag = freadshow(0, fid,1,'*uint32');
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_UNDEFINED = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00000000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_POLYGRAPHIC = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00010000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_TRIGGER = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00020000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_CORTICALGRID = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00040000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_INTRACRANIAL = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00080000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_SCALPELECTRODE = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00100000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_MAGNETOMETER = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00200000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_AXIAL_GRADIOMETER = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00400000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_PLANAR_GRADIOMETER = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('01000000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_MEGREFERENCE = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00800000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_NKC_REFERENCE = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('02000000')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANTYPE_CHANSTATE_BAD = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00000001')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANSTATE_REFERENCE = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00000002')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANSTATE_INTERPOLRECORDED = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00000004')),'uint32'));
      BCAL_struct.channel_states(channel_n).BSA_CHANSTATE_INVISIBLE = logical(bitand(BCAL_struct.channel_states(channel_n).flag,uint32(hex2dec('00001000')),'uint32'));
    case 'CHCO'
      % Channel coordinates in mm
      n_channels = current_length / 4 / 9; % Divide by 4 for *single and by 9 for number of elements per channel
      BCAL_struct.coord_data = zeros(n_channels,9,'single');
      for channel_n = 1:n_channels
        BCAL_struct.coord_data(channel_n,:) = freadshow(0, fid,9,'*single');
      end
      % More processing done later to obtain actual coordinates
    case 'CHNR'
      % Total number of channels
      BCAL_struct.n_channels = double(freadshow(0, fid,1,'*uint16'));
    case 'CHLA'
      % Channel label of a channel with the specified index
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.channel_labels{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHET'
      % Electrode thickness
      BCAL_struct.electrode_thickness = freadshow(0, fid,1,'*single');
    case 'CHSI'
      % Spline interpolation smoothing constant
      BCAL_struct.spline_smoothing_constant = freadshow(0, fid,1,'*single');
    case 'CHLS'
      % Least significant bits of data
      BCAL_struct.lsbs = double(freadshow(0, fid,current_length/4,'*single'));
      % Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead. %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CHSF'
      % Sampling frequency
      BCAL_struct.s_rates = freadshow(0, fid,current_length/4,'*double');
    case 'HCMM'
      % Head center in mm
      BCAL_struct.head_center.x = freadshow(0, fid,1,'*single');
      BCAL_struct.head_center.y = freadshow(0, fid,1,'*single');
      BCAL_struct.head_center.z = freadshow(0, fid,1,'*single');
    case 'HRMM'
      % Head radius in mm
      BCAL_struct.head_radius = freadshow(0, fid,1,'*single');
    case 'FIDC'
      % Fiducial coordinates in mm
      BCAL_struct.fiducial.nasion.x = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.nasion.y = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.nasion.z = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.left_preauricular.x = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.left_preauricular.y = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.left_preauricular.z = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.right_preauricular.x = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.right_preauricular.y = freadshow(0, fid,1,'*single');
      BCAL_struct.fiducial.right_preauricular.z = freadshow(0, fid,1,'*single');
    case 'HSPN'
      % Total number of head surface points
      BCAL_struct.n_addn_surf_pnts = double(freadshow(0, fid,1,'*int16'));
    case 'HSPC'
      % Head surface point coordinates
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.head_surface_points{channel_n}.x = freadshow(0, fid,1,'*single');
      BCAL_struct.head_surface_points{channel_n}.y = freadshow(0, fid,1,'*single');
      BCAL_struct.head_surface_points{channel_n}.z = freadshow(0, fid,1,'*single');
    case 'HSPD'
      % Head surface point labels
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.head_surface_points{channel_n}.label = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHCU'
      % Channel units
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.channel_units{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
    case 'CHFI'
      % Filter information
      fseek(fid,current_length,'cof');
      fprintf('Skipping\t%s\n', current_tag);%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      % SKIP FOR NOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      % Somewhat complicated structure, no test data %%%%%%%%%%%%%%%%%% TODO
    case 'CHNU'
      % Channel numbers
      BCAL_struct.channel_ns = freadshow(0, fid,current_length/4,'*int32');
    case 'CHCM'
      % Channel comments
      channel_n = double(freadshow(0, fid,1,'*uint16'))+1; % Plus 1 because index starts at 0
      BCAL_struct.channel_comments{channel_n} = read_chars(fid,current_length-2); % Subtract 2 from offet for channel_n
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
if((BCAL_struct.offsets(end)+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from channel block. Should have read %d bytes', ...
    (ftell(fid)-BCAL_struct.offsets(end))-expected_length,ftell(fid)-BCAL_struct.offsets(end),expected_length);
end


%% EVENT BLOCK FUNCTIONS

function BEVT_struct = read_BEVT(fid, file_length, BEVT_struct, BEVT_offset)
%% Read event block
% BEVT_offset [scalar] - offset of current event block start
% BEVT_struct [structure] - Existing or blank BEVT structure
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
      BEVT_struct.n_events = double(freadshow(1, fid,1,'*uint32'));
    case 'VERS'
      BEVT_struct.version = double(freadshow(1, fid,1,'*uint32'));
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
for event_n = 1:BEVT_struct.n_events
  [current_tag,current_length] = read_tag_offset_pair(fid);
  BEVT_struct.events{event_n} = read_event_tag(fid,current_tag,current_length);
end

% Check that expected amout of file was read
expected_length = double(event_block_length) + 8; % 8 for tag and offset
if((BEVT_offset+expected_length) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from event block. Should have read %d bytes', ...
    (ftell(fid)-BEVT_offset)-expected_length,ftell(fid)-BEVT_offset,expected_length);
end

function event_obj = read_event_tag(fid,event_tag,event_length)
% Read an event into a structure

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Create the event object
event_obj.evttype = event_tag;

switch event_tag
  case 'BASE'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Base event tag
    event_obj = read_event_tag_base(fid,ftell(fid),event_length,event_obj);
  case 'COMM'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Comment event tag
    event_obj = read_event_tag_comm(fid,ftell(fid),event_length,event_obj);
  case 'MARK'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Marker event tag
    event_obj = read_event_tag_mark(fid,ftell(fid),event_length,event_obj);
  case 'GENE'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Generic event tag
    event_obj = read_event_tag_gene(fid,ftell(fid),event_length,event_obj);
  case 'SEGM'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Segment event tag
    event_obj = read_event_tag_segm(fid,ftell(fid),event_length,event_obj);
  case 'ASGM'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Average segment start event tag
    event_obj = read_event_tag_asgm(fid,ftell(fid),event_length,event_obj);
  case 'MPS'
    % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPS]',event_length);
    end
  case 'MPSC'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Classified multiple pattern search event tag
    %   used by BESA internally
    if(fseek(fid,event_length,'cof') == -1)
      fclose(fid);
      error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [LIST:MPSC]',event_length);
    end
  case 'PATT'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Pattern event tag
    event_obj = read_event_tag_patt(fid,ftell(fid),event_length,event_obj);
  case 'TRIG'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Trigger event tag
    event_obj = read_event_tag_trig(fid,ftell(fid),event_length,event_obj);
  case 'PAIR'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Paired event tag
    event_obj = read_event_tag_pair(fid,ftell(fid),event_length,event_obj);
  case 'ARTI'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Artifact event tag
    event_obj = read_event_tag_arti(fid,ftell(fid),event_length,event_obj);
  case 'EPOC'
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
    % Epoch event tag
    event_obj = read_event_tag_epoc(fid,ftell(fid),event_length,event_obj);
  case 'IMP'
    % THREE CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    if verbose; fprintf('%s %d\n',event_tag,event_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (base_offset+base_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SAMP'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Sample index (zero based)
      event_obj.sample_n = freadshow(0, fid,1,'*int64');
    case 'TIME'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Event time
      event_obj.time.year = freadshow(0, fid,1,'*int16');
      event_obj.time.month = freadshow(0, fid,1,'*int16');
      event_obj.time.dayOfWeek = freadshow(0, fid,1,'*int16');
      event_obj.time.day = freadshow(0, fid,1,'*int16');
      event_obj.time.hour = freadshow(0, fid,1,'*int16');
      event_obj.time.minute = freadshow(0, fid,1,'*int16');
      event_obj.time.second = freadshow(0, fid,1,'*int16');
      event_obj.time.milliseconds = freadshow(0, fid,1,'*int16');
      event_obj.time.microseconds = freadshow(0, fid,1,'*single');
      event_obj.time.stateFlag = freadshow(0, fid,1,'*uint64');
    case 'SIDX'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Segment index (zero based)
      event_obj.segment_index = freadshow(0, fid,1,'*int32');
    case 'CODE'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Event code (zero based)
      % This value is used by events of type Pattern (PATT) and Trigger (TRIG)
      %  to store the pattern number and the trigger code. Note that the number/code minus 1 is stored.
      %  Additionally, events of type Artifact (ARTI) and Epoch (EPOC) use the code value
      %  internally due to historical reasons. Other event types may use the code value to
      %  store additional information.
      event_obj.code = freadshow(0, fid,1,'*int32');
    case 'EVID'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Internal BESA event ID (zero based)
      event_obj.besa_event_id = freadshow(0, fid,1,'*int32');
    case 'STAT'
      if verbose; fprintf('BASE:%s %d\n',current_tag,current_length); end;
      % Event state
      event_obj.state.value = freadshow(0, fid,1,'*uint32');
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (comm_offset+comm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'TEXT'
      if verbose; fprintf('COMM:%s %d\n',current_tag,current_length); end;
      % Event text
      event_obj.text = read_chars(fid,current_length);
    case 'BASE'
      if verbose; fprintf('COMM:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (mark_offset+mark_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      if verbose; fprintf('MARK:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (gene_offset+gene_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'COMM'
      if verbose; fprintf('GENE:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (segm_offset+segm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'SBEG'
      if verbose; fprintf('SEGM:%s %d\n',current_tag,current_length); end;
      % Segment start time
      event_obj.segment_start.year = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.month = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.dayOfWeek = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.day = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.hour = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.minute = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.second = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.milliseconds = freadshow(0, fid,1,'*int16');
      event_obj.segment_start.microseconds = freadshow(0, fid,1,'*single');
      event_obj.segment_start.stateFlag = freadshow(0, fid,1,'*uint64');
    case 'DAYT'
      if verbose; fprintf('SEGM:%s %d\n',current_tag,current_length); end;
      % Day time of segment start in microseconds
      event_obj.segment_start.dayt = freadshow(0, fid,1,'*double');
    case 'INTE'
      if verbose; fprintf('SEGM:%s %d\n',current_tag,current_length); end;
      % Sampling interval in microseconds
      event_obj.segment_start.sampling_interval = freadshow(0, fid,1,'*double');
    case 'COMM'
      if verbose; fprintf('SEGM:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (asgm_offset+asgm_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'STIM'
      if verbose; fprintf('ASGM:%s %d\n',current_tag,current_length); end;
      % Prestimulus baseline interval in microseconds
      event_obj.baseline_interval = freadshow(0, fid,1,'*double');
    case 'AVRS'
      if verbose; fprintf('ASGM:%s %d\n',current_tag,current_length); end;
      % Number of averages
      event_obj.n_averages = freadshow(0, fid,1,'*int32');
    case 'COMM'
      if verbose; fprintf('ASGM:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (patt_offset+patt_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'BASE'
      if verbose; fprintf('PATT:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (trig_offset+trig_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'COMM'
      if verbose; fprintf('TRIG:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

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
      if verbose; fprintf('PAIR:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (arti_offset+arti_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      if verbose; fprintf('ARTI:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (epoc_offset+epoc_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'PAIR'
      if verbose; fprintf('EPOC:%s %d\n',current_tag,current_length); end;
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

% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Loop through all tags in data section
while ~feof(fid) && ftell(fid) < (imp_offset+imp_length)
  [current_tag,current_length] = read_tag_offset_pair(fid);
  switch current_tag
    case 'FORM'
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
      % Indicates the format used to store impedance values in VAL
      %  Set to 0 if the impedance status (valid/invalid) is stored
      %  Set to 1 if impedance values (in kOhm) are stored
      event_obj.impedance.format = freadshow(0, fid,1,'*int32');
    case 'NR'
      % TWO CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
      % Number of channels for which impedance information is stored.
      %   (Number of elements stored in TYPE, LABL and VAL)
      event_obj.impedance.n_channels = freadshow(0, fid,1,'*uint32');
    case 'TYPE'
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
      % Channel types
      % The flags used for channel type description are the same as used for the CHTS data elements (specified in chapter 2.3).
      %   Note: Channel type and channel label are used for identification of channel
      %   for which impedance information is set. (Compare to data elements CHTS and CHLA, specified in chapter 2.3).
      event_obj.impedance.types = freadshow(0, fid,event_obj.impedance.n_channels,'*uint32'); % Assumes that n_channels has been set
    case 'LABL'
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
      % Channel labels
      %   Note: Channel type and channel label are used for identification of channel for which impedance information is set.
      %   (Compare to data elements CHTS and CHLA as specified in chapter 2.3).
      event_obj.impedance.labels = read_chars(fid,current_length);
    case 'VAL'
      % TWO CHARACTERS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
      % Impedance values
      % Depending on format set in FORM either an impedance STATUS (ok/not ok) or an impedance VALUE (in kOhm) is stored
      %   A value of -1 means that the impedance is not set or invalid
      event_obj.impedance.values = freadshow(0, fid,event_obj.impedance.n_channels,'*single'); % Assumes that n_channels has been set
    case 'BASE'
      if verbose; fprintf('IMP:%s %d\n',current_tag,current_length); end;
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

%% DATA BLOCK FUNCTIONS

function block_data = read_BDAT(fid, file_length, bdat_offset, n_channels, lsbs)
%% Read data block
%
% lsbs - [1 x n_channels array] - int data is multiplied by this for scaling

% CHECK FOR LSB LESS THAN ZERO %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
% Please note that zero or negative LSB values are not allowed. If a non-positive value is found in the array, a value of "1.f" is used instead. %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% REMOVE THIS AND ALL FPRINTF STATEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
verbose = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO

% Constants for data type
CONST_FLOAT = 0;
CONST_INT16 = 1;
CONST_COMPRESSED = 1;
CONST_UNCOMPRESSED = 0;
% Constants for prefix byte
CONST_NOCOMPRESSION = 0;
CONST_FIRSTSCHEME = 3;
CONST_SECONDSCHEME = 4;
CONST_THIRDSCHEME = 5;
CONST_NOCOMPRESSION_FIRST2INT = 6;
CONST_FIRSTSCHEME_FIRST2INT = 7;
CONST_NOCOMPRESSION_ALLINT = 8;
CONST_ZLIB_DD = 9;
CONST_ZLIB_FIRSTSCHEME = 13;
CONST_ZLIB_SECONDSCHEME = 14;
CONST_ZLIB_THIRDSCHEME = 15;
CONST_ZLIB_FIRSTSCHEME_FIRST2INT = 17;
CONST_ZLIB_SECONDSCHEME_FIRST2INT = 18;
CONST_ZLIB_THIRDSCHEME_FIRST2INT = 19;
CONST_ZLIB_DD_ALLINT = 29;

% Skip to start of BDAT section
if(fseek(fid,double(bdat_offset),'bof') == -1) % double() because Windows can't seek to uint64
  fclose(fid);
  error('ReadBesaMatlab:ErrorFseek','fseek to %d failed [BDAT]',bdat_offset);
end

% Read BDAT tag and offset
[~,ofst_BDAT] = read_tag_offset_pair(fid,'BDAT');

% Check that file is not shorter than expected
if(file_length < (ftell(fid)+ofst_BDAT))
  expected_length = ftell(fid)+ofst_BDAT;
  fclose(fid);
  error('ReadBesaMatlab:ErrorFileTooShortForDataBlock','Data block expected file at least %d bytes long but file is %d bytes long',expected_length,file_length);
end

% Determine type of data in this block
read_tag_offset_pair(fid,'DATT');
flag_BDAT = freadshow(1, fid,1,'*uint32');
if bitand(flag_BDAT,uint32(1),'uint32')
  data_type = CONST_INT16;
else
  data_type = CONST_FLOAT;
end
if bitand(flag_BDAT,uint32(hex2dec('0010')),'uint32')
  data_comp = CONST_COMPRESSED;
else
  data_comp = CONST_UNCOMPRESSED;
end

% Determine number of samples in this block
read_tag_offset_pair(fid,'DATS');
n_samples = freadshow(1, fid,1,'*uint32');

% Read DATA tag and offset
[~,data_block_length] = read_tag_offset_pair(fid,'DATA');
data_block_offset = double(ftell(fid));

% Read data
block_data = zeros(n_samples,n_channels,'double');
switch num2str([data_type data_comp])
  
  case num2str([CONST_INT16 CONST_UNCOMPRESSED])
    % Read int16s, reshape to [n_samples x n_channels], multiply each channel by LSB
    block_data = bsxfun(@times,lsbs', ...
      double(reshape(freadshow(1, fid,n_samples*n_channels,'*int16'),[n_samples,n_channels])));
    
  case num2str([CONST_FLOAT CONST_UNCOMPRESSED])
    % Read singles, reshape to [n_samples x n_channels]
    block_data = double(reshape(freadshow(1, fid,n_samples*n_channels,'*single'),[n_samples,n_channels]));
    
  case {num2str([CONST_FLOAT CONST_COMPRESSED]),num2str([CONST_INT16 CONST_COMPRESSED])}
    % Compressed data
    for channel_n = 1:n_channels
      if verbose; fprintf('PREFIX VALUE:\n'); end; %%%%%%%%%%%%%%%%%% TODO
      prefix_val = freadshow(1, fid,1,'*uint8');
      switch prefix_val
        case CONST_NOCOMPRESSION
          if verbose; fprintf('No Compression\n'); end;
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int16. Rest are int16.
          block_data(:,channel_n) = double(freadshow(1, fid,n_samples,'*int16'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_FIRSTSCHEME
          if verbose; fprintf('First scheme\n'); end;
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          first_2_vals = freadshow(1, fid,2,'*int16');
          block_data(:,channel_n) = double(decode_firstscheme(fid,freadshow(1, fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_SECONDSCHEME
          if verbose; fprintf('Second scheme, first two values:\n'); end;
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, second scheme)
          first_2_vals = freadshow(1, fid,2,'*int16');
          block_data(:,channel_n) = double(decode_secondscheme(fid,freadshow(1, fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_THIRDSCHEME
          if verbose; fprintf('Third scheme\n'); end;
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, third scheme)
          first_2_vals = freadshow(1, fid,2,'*int16');
          block_data(:,channel_n) = double(decode_thirdscheme(fid,freadshow(1, fid,data_block_length-4,'*uint8'), n_samples, first_2_vals));
        case CONST_NOCOMPRESSION_FIRST2INT
          if verbose; fprintf('No comp first 2 int\n'); end;
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int16
          block_data(1:2,channel_n) = double(freadshow(1, fid,2,'*int32'));
          block_data(3:end,channel_n) = double(freadshow(1, fid,n_samples-2,'*int16'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_FIRSTSCHEME_FIRST2INT
          if verbose; fprintf('First scheme first 2 int\n'); end;
          % No zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, first scheme)
          first_2_vals = freadshow(1, fid,2,'*int32');
          block_data(:,channel_n) = double(decode_firstscheme(fid,freadshow(1, fid,data_block_length-8,'*uint8'), n_samples, first_2_vals));
        case CONST_NOCOMPRESSION_ALLINT
          if verbose; fprintf('No compression all int\n'); end;
          % No zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int32
          block_data(:,channel_n) = double(freadshow(1, fid,n_samples,'*int32'));
          % Integrate twice
          block_data(2:end,channel_n) = cumsum(block_data(2:end,channel_n),1);
          block_data(:,channel_n) = cumsum(block_data(:,channel_n),1);
        case CONST_ZLIB_DD
          if verbose; fprintf('Zlib dd\n'); end;
          % Yes zlib. No pre-compression. Yes double difference
          % First two elements are int16. Rest are int16
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = typecast(buffer_data,'int16');
        case CONST_ZLIB_FIRSTSCHEME
          if verbose; fprintf('Zlib first scheme\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_firstscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_SECONDSCHEME
          if verbose; fprintf('Zlib second scheme\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, second scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_secondscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_THIRDSCHEME
          if verbose; fprintf('Zlib third scheme\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, third scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_thirdscheme(fid,buffer_data(5:end), n_samples, typecast(buffer_data(1:4),'int16')));
        case CONST_ZLIB_FIRSTSCHEME_FIRST2INT
          if verbose; fprintf('Zlib first scheme first 2 int\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int16. Rest are int8 (pre-compressed, first scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_firstscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_SECONDSCHEME_FIRST2INT
          if verbose; fprintf('Zlib second scheme first 2 int\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, second scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_secondscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_THIRDSCHEME_FIRST2INT
          if verbose; fprintf('Zlib third scheme first 2 int\n'); end;
          % Yes zlib. Yes pre-compression. Yes double difference
          % First two elements are int32. Rest are int8 (pre-compressed, third scheme)
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*uint8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = double(decode_thirdscheme(fid,buffer_data(9:end), n_samples, typecast(buffer_data(1:8),'int32')));
        case CONST_ZLIB_DD_ALLINT
          if verbose; fprintf('Zlib dd all int\n'); end;
          % Yes zlib. No pre-compression. Yes double difference
          % First two elements are int32. Rest are int32
          buffer_len = freadshow(1, fid,1,'*uint32');
          buffer_data = freadshow(1, fid,buffer_len,'*int8');
          buffer_data = typecast(dunzip(buffer_data),'uint8')';
          if verbose; fprintf('\tUnzipped to %d bytes\n',numel(buffer_data)); end;
          block_data(:,channel_n) = typecast(buffer_data,'int32'); 
        otherwise
          current_loc = ftell(fid);
          fclose(fid);
          error('ReadBesaMatlab:ErrorBDATReadPrefixValueUnknownScheme','Unknown scheme  CH:%d  prefix_val:%d  File offset:%d',channel_n,prefix_val,current_loc);
      end
    end
    
    if(strcmp(num2str([data_type data_comp]),num2str([CONST_INT16 CONST_COMPRESSED])))
      % Multiply int16 data by lsbs
      block_data = bsxfun(@times,lsbs',block_data);
      % CHECK FOR LSB LESS THAN ZERO %%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    end
end

% Check that expected amout of data was read
if((data_block_offset+double(data_block_length)) ~= ftell(fid))
  warning('ReadBesaMatlab:WarningDidNotReadExactBlockLength','%d bytes off. Read %d bytes from data block. Should have read %d bytes', ...
    (ftell(fid)-data_block_offset)-double(data_block_length),ftell(fid)-data_block_offset,double(data_block_length));
end

function outbuffer = decode_firstscheme(fid, inbuffer, n_samples, first2vals)
% Read data in first scheme

CONST_MESHGRID_VALS_1 = -7:7;
CONST_AB_INT32_RANGE = 241:-1:236; % Reverse order. This is needed to determine n_vals
CONST_AB_INT16_RANGE = 247:-1:242;
CONST_AB_INT8_RANGE  = 254:-1:248;

max_lut_val = numel(CONST_MESHGRID_VALS_1)^2-1; % Any buffer value greater than this is an announcing byte

% Use persistent variable so lookup table does not need to be recomputed each time
persistent firstscheme_lookuptable;
if isempty(firstscheme_lookuptable)
  % Create the lookup grid from -7 to 7 in x and y
  [firstscheme_lookuptable(:,:,1),firstscheme_lookuptable(:,:,2)] = meshgrid(CONST_MESHGRID_VALS_1,CONST_MESHGRID_VALS_1);
  % Reshape the lookup grid to be [1:225 x 1:2]
  firstscheme_lookuptable = reshape(firstscheme_lookuptable,[numel(CONST_MESHGRID_VALS_1)^2 2]);
end

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  %   Transpose and then use linear indexing in the output to put all
  %   elements into a 1-d array
  try
    outbuffer((last_outbuffer_idx+1):end) = firstscheme_lookuptable(inbuffer+1); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(firstscheme_lookuptable(inbuffer+1));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+2*(ab_idx-last_ab_idx-1))) = ...
      firstscheme_lookuptable(inbuffer((last_ab_idx+1):(ab_idx-1))+1,:); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+2*(ab_idx-last_ab_idx-1))));
      received_samples = numel(firstscheme_lookuptable(inbuffer((last_ab_idx+1):(ab_idx-1))+1,:));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, middle of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
  last_outbuffer_idx = (last_outbuffer_idx+2*(ab_idx-last_ab_idx-1));
  
  if(any(CONST_AB_INT32_RANGE == inbuffer(ab_idx)))
    % AB indicates int32
    n_vals = find(CONST_AB_INT32_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*4; % x4 for int32
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int32');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an alowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      firstscheme_lookuptable(inbuffer((last_ab_idx+1):end)+1,:); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(firstscheme_lookuptable(inbuffer((last_ab_idx+1):end)+1,:));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [first scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Integrate twice
outbuffer(2:end) = cumsum(outbuffer(2:end));
outbuffer = cumsum(outbuffer);

function outbuffer = decode_secondscheme(fid, inbuffer, n_samples, first2vals)
% Decode second scheme

CONST_MESHGRID_VALS_2A = -2:2;
CONST_MESHGRID_VALS_2B = -5:5;
CONST_AB_INT16_RANGE = 249:-1:246; % Reverse order. This is needed to determine n_vals
CONST_AB_INT8_RANGE  = 254:-1:250;
meshgrid_vals.A = CONST_MESHGRID_VALS_2A;
meshgrid_vals.B = CONST_MESHGRID_VALS_2B;

max_lut_val = numel(CONST_MESHGRID_VALS_2A)^3 + numel(CONST_MESHGRID_VALS_2B)^2 - 1; % Any buffer value greater than this is an announcing byte

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  try
    outbuffer((last_outbuffer_idx+1):end) = secondscheme_lookup(inbuffer+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(secondscheme_lookup(inbuffer+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [second scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  % No error checking, because we don't know how long it should be
  decoded_buffer = secondscheme_lookup(inbuffer((last_ab_idx+1):(ab_idx-1))+1,meshgrid_vals); % Plus 1 because indices start at 0
  outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+numel(decoded_buffer))) = ...
    decoded_buffer;
  last_outbuffer_idx = (last_outbuffer_idx+numel(decoded_buffer));
  clear decoded_buffer;
  
  if(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an allowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range [second scheme]: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      secondscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(secondscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [second scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

function output = secondscheme_lookup(input, meshgrid_vals)
% Lookup table for second scheme

% Use persistent variable so lookup table does not need to be recomputed each time
persistent secondscheme_lookuptable;
if isempty(secondscheme_lookuptable)
  
  % Create the lookup grid from -2 to 2 in x, y, z
  [secondscheme_lookuptable_a(:,:,:,1),secondscheme_lookuptable_a(:,:,:,2),secondscheme_lookuptable_a(:,:,:,3)] = ...
    meshgrid(meshgrid_vals.A,meshgrid_vals.A,meshgrid_vals.A);
  % Reshape the lookup grid to be [1:125 x 1:3]
  secondscheme_lookuptable_a = reshape(secondscheme_lookuptable_a,[numel(meshgrid_vals.A)^3 3]);
  % Correct order of x,y,z
  secondscheme_lookuptable_a(:,[1 2 3]) = secondscheme_lookuptable_a(:,[3 1 2]);
  
  % Create the lookup grid from -5 to 5 in x and y
  [secondscheme_lookuptable_b(:,:,1),secondscheme_lookuptable_b(:,:,2)] = meshgrid(meshgrid_vals.B,meshgrid_vals.B);
  % Reshape the lookup grid to be [1:121 x 1:2]
  secondscheme_lookuptable_b = reshape(secondscheme_lookuptable_b,[numel(meshgrid_vals.B)^2 2]);
  
  % Put the lookup tables together in a cell array (because of different sized cells)
  secondscheme_lookuptable = num2cell(secondscheme_lookuptable_a,2);
  secondscheme_lookuptable = [secondscheme_lookuptable; num2cell(secondscheme_lookuptable_b,2)];
  
  clear secondscheme_lookuptable_a;
  clear secondscheme_lookuptable_b;
end

output_cell = secondscheme_lookuptable(input);
output = [output_cell{:}];

function outbuffer = decode_thirdscheme(fid, inbuffer, n_samples, first2vals)
% Decode third scheme

CONST_MESHGRID_VALS_3A = -1:1;
CONST_MESHGRID_VALS_3B = -6:6;
CONST_AB_INT16_RANGE = 251:-1:250; % Reverse order. This is needed to determine n_vals
CONST_AB_INT8_RANGE  = 254:-1:252;
meshgrid_vals.A = CONST_MESHGRID_VALS_3A;
meshgrid_vals.B = CONST_MESHGRID_VALS_3B;

max_lut_val = numel(CONST_MESHGRID_VALS_3A)^4 + numel(CONST_MESHGRID_VALS_3B)^2 - 1; % Any buffer value greater than this is an announcing byte

% Initialize outbuffer
outbuffer = zeros(n_samples,1,'int32');

% Fill in the first two values
outbuffer(1:2) = first2vals;

% Find first announcing byte (AB) (value outside of LUT)
ab_idx = find(inbuffer>max_lut_val,1,'first');

last_outbuffer_idx = 2; % first2vals
if isempty(ab_idx)
  % No ABs, just use lookup table for whole inbuffer
  % Get the output from the lookup table
  try
    outbuffer((last_outbuffer_idx+1):end) = thirdscheme_lookup(inbuffer+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(thirdscheme_lookup(inbuffer+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [third scheme, no ABs]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

% Loop until out of announcing bytes
possible_abs = inbuffer > max_lut_val;
last_ab_idx = 0;
while ~isempty(ab_idx)
  
  % Fill outbuffer using LUT with all values between the last set of non-encodable values 
  %   and the current set of non-encodable values,
  %   starting at the last filled outbuffer index.
  % No error checking, because we don't know how long it should be
  decoded_buffer = thirdscheme_lookup(inbuffer((last_ab_idx+1):(ab_idx-1))+1,meshgrid_vals); % Plus 1 because indices start at 0
  outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+numel(decoded_buffer))) = ...
    decoded_buffer;
  last_outbuffer_idx = (last_outbuffer_idx+numel(decoded_buffer));
  clear decoded_buffer;
  
  if(any(CONST_AB_INT16_RANGE == inbuffer(ab_idx)))
    % AB indicates int16
    n_vals = find(CONST_AB_INT16_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals*2; % x2 for int16
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int16');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  elseif(any(CONST_AB_INT8_RANGE == inbuffer(ab_idx)))
    % AB indicates int8
    n_vals = find(CONST_AB_INT8_RANGE==inbuffer(ab_idx),1);
    n_skip = n_vals; % x1 for int8
    % Fill outbuffer with n_vals
    outbuffer((last_outbuffer_idx+1):(last_outbuffer_idx+n_vals)) = typecast(inbuffer((ab_idx+1):(ab_idx+n_skip)),'int8');
    last_outbuffer_idx = last_outbuffer_idx+n_vals;
    last_ab_idx = ab_idx+n_skip;
  else
    % not an allowed announcing byte value
    fclose(fid);
    error('ReadBesaMatlab:ErrorABOutOfRange','Announcing byte out of range [third scheme]: %d',inbuffer(ab_idx));
  end
  
  % Go to next AB
  ab_idx = last_ab_idx + find(possible_abs((last_ab_idx+1):end),1,'first'); % Note: X+[]=[]
  
end

if(last_ab_idx<numel(inbuffer))
  % Fill outbuffer using LUT with all values after the last set of non-encodable values
  %   starting at the last filled outbuffer index.
  try
    outbuffer((last_outbuffer_idx+1):end) = ...
      thirdscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals); % Plus 1 because indices start at 0
  catch ME
    if(strcmp(ME.identifier,'MATLAB:subsassignnumelmismatch'))
      expected_samples = numel(outbuffer((last_outbuffer_idx+1):end));
      received_samples = numel(thirdscheme_lookup(inbuffer((last_ab_idx+1):end)+1,meshgrid_vals));
      fclose(fid);
      error('ReadBesaMatlab:ErrorUnexpectedNSamplesFromPreCompression','Expected %d samples, but got %d samples. [third scheme, end of buffer]', ...
        expected_samples,received_samples);
    else
      rethrow(ME);
    end
  end
end

function output = thirdscheme_lookup(input, meshgrid_vals)
% Lookup table for third scheme

% Use persistent variable so lookup table does not need to be recomputed each time
persistent thirdscheme_lookuptable;
if isempty(thirdscheme_lookuptable)
  
  % Create the lookup grid from -1 to 1 in x, y, z, c
  [thirdscheme_lookuptable_a(:,:,:,:,1),thirdscheme_lookuptable_a(:,:,:,:,2),thirdscheme_lookuptable_a(:,:,:,:,3),thirdscheme_lookuptable_a(:,:,:,:,4)] = ...
    ndgrid(meshgrid_vals.A);
  % Reshape the lookup grid to be [1:81 x 1:4]
  thirdscheme_lookuptable_a = reshape(thirdscheme_lookuptable_a,[numel(meshgrid_vals.A)^4 4]);
  % Correct order of x,y,z,c
  thirdscheme_lookuptable_a(:,[1 2 3 4]) = thirdscheme_lookuptable_a(:,[4 3 2 1]);
  
  % Create the lookup grid from -6 to 6 in x and y
  [thirdscheme_lookuptable_b(:,:,1),thirdscheme_lookuptable_b(:,:,2)] = meshgrid(meshgrid_vals.B,meshgrid_vals.B);
  % Reshape the lookup grid to be [1:169 x 1:2]
  thirdscheme_lookuptable_b = reshape(thirdscheme_lookuptable_b,[numel(meshgrid_vals.B)^2 2]);
  
  % Put the lookup tables together in a cell array (because of different sized cells)
  thirdscheme_lookuptable = num2cell(thirdscheme_lookuptable_a,2);
  thirdscheme_lookuptable = [thirdscheme_lookuptable; num2cell(thirdscheme_lookuptable_b,2)];
  
  clear thirdscheme_lookuptable_a;
  clear thirdscheme_lookuptable_b;
end

output_cell = thirdscheme_lookuptable(input);
output = [output_cell{:}];


%% HELPER FUNCTIONS

function M = dunzip(Z)
% DUNZIP - decompress gzipped stream of bytes
% FORMAT M = dzip(Z)
% Z  -  compressed variable to decompress (uint8 vector)
% M  -  decompressed output
% 
% See also DZIP

% Carefully tested, but no warranty; use at your own risk.
% Michael Kleder, Nov 2005
% Modified by Guillaume Flandin, May 2008

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
a   = java.io.ByteArrayInputStream(Z);
b   = java.util.zip.InflaterInputStream(a);
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
c   = java.io.ByteArrayOutputStream;
isc.copyStream(b,c);
M   = c.toByteArray;



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
fprintf([sprintf('%.15d',ftell(fid)) ': ']); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
fprintf('TAG: [%s] ',out_tag); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
fprintf('OFFSET: %d\n',out_offset); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO


function outchars = read_chars(fid,n_chars)
% Read n_chars characters from file at current position
%   Replace null character (aka char(0) or '\x0') with ''
%   Note transpose after fread
outchars = regexprep(freadshow(1, fid,n_chars,'*char')','\x0','');

function dataout = freadshow(printoptions,fid,n_items,vartype)
% REMOVE THIS AND REPLACE ALL CALLS WITH fread() %%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
% fread() and print output (for testing)
% printoptions:
%  0 - Don't print
%  1 - Print and don't pause
%  2 - Print and pause
pauseon = printoptions - 1;
if printoptions
  fprintf([sprintf('% 15d',ftell(fid)) ': ']);
end
dataout = fread(fid,n_items,vartype);
if printoptions
  if(n_items<10)
    if(ischar(dataout))
      try
        fprintf([dataout '\n']);
      catch
        fprintf([dataout' '\n']);
      end
    else
      for i = 1:n_items
        fprintf([num2str(dataout(i)) ' ']);
      end
      fprintf('\n');
    end
  else
    fprintf('**** %d %s read ****\n',n_items,vartype);
  end
  if pauseon
    pause;
  end
end












