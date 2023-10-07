function N = spm_dicom_metadata(N,hdr)
% Export image metadata as side-car JSON file
% FORMAT N = spm_dicom_metadata(N,hdr)
% N(input)  - nifti object
% hdr       - a single header from spm_dicom_headers
% N(output) - unchanged nifti object (for potential future use)
%
% This function creates JSON-encoded metadata during DICOM to NIfTI
% conversion, including all acquisition parameters, and saves them as a
% JSON side-car file.
%
% See also: spm_dicom_convert
%__________________________________________________________________________
% Copyright (C) 2017-2018 Wellcome Trust Centre for Neuroimaging

% Evelyne Balteau, Cyclotron Research Centre, University of Liege, Belgium
% $Id: spm_dicom_metadata.m 7482 2018-11-12 12:18:08Z guillaume $


%-Provenance and general description of the image (hMRI toolbox)
dicom_convert_version = {sprintf('%s %s', spm_check_version, version),sprintf('spm_dicom_convert.m - %s', spm('Version'))};
metadata.history.procstep = struct('descrip','dicom to nifti import', 'version', {dicom_convert_version}, 'procpar', []);
metadata.history.input(1) = struct('filename','AnonymousFileName', 'history',[]);
if isfield(hdr,'ImageType')
    metadata.history.output = struct('imtype',hdr.ImageType, 'units','a.u.');
else
    metadata.history.output = struct('imtype','Unprocessed MR image', 'units','a.u.');
end

%-Acquisition parameters with complete DICOM header (hMRI toolbox)
hdr = reformat_spm_dicom_header(hdr);
hdr = anonymise_metadata(hdr); % default is basic anonymisation
metadata.acqpar = hdr;

%-Save metadata as side-car JSON file
spm_jsonwrite(spm_file(N.dat.fname,'ext','json'), metadata, struct('indent','\t'));


%==========================================================================
function hdr = reformat_spm_dicom_header(hdr)
% To tidy up and rearrange CSA fields in the header, including formatting
% the ASCII part into a proper MATLAB structure (Note: this is specific to
% Siemens DICOM format but could be extended to other cases where needed)
%
% FORMAT hdr = reformat_spm_dicom_header(hdr)
% where hdr is a dicom header structure as output by spm_dicom_headers.

% list of CSA-like fields:
CSAlist = {'CSAImageHeaderInfo', 'CSASeriesHeaderInfo','CSANonImageHeaderInfoVA','CSAMiscProtocolHeaderInfoVA','CSANonImageHeaderInfoVB','CSAMiscProtocolHeaderInfoVB'};
% NB: might be necessary to add cases 'Private_0029_1110' and
% 'Private_0029_1210' for spectroscopic data (see spm_dicom_essentials.m)

for ccsa = 1:length(CSAlist)
    if isfield(hdr,CSAlist{ccsa})
        tdyhdr = tidy_CSA(hdr.(CSAlist{ccsa}));
        if isfield(tdyhdr,'MrPhoenixProtocol')
            % works for Siemens VB, VD & VE DICOM format:
            tdyhdr.MrPhoenixProtocol = read_ASCII(tdyhdr.MrPhoenixProtocol);
        elseif isfield(tdyhdr,'MrProtocol')
            % works for Siemens VA DICOM format:
            tdyhdr.MrProtocol = read_ASCII(tdyhdr.MrProtocol);
        end
        hdr.(CSAlist{ccsa}) = tdyhdr;
    end
end


%==========================================================================
function hdrout = anonymise_metadata(hdr, opts)
% FORMAT hdrout = anonymise_metadata(hdr, opts)
% hdr is a DICOM header read with spm_dicom_header
% opts are anonymisation options:
%       opts.anonym = 'none': no anonymisation, all patient data kept in
%                       the metadata.
%                     'full': no patient information is kept at all
%                     'basic': patient ID (presumably not containing his
%                       name), age (years at the time of the data
%                       acquisition), sex, size and weight are kept,
%                       patient name, date of birth and DICOM filename
%                       (often containing the patient name) are removed.
%
% !!! IMPORTANT WARNING: EFFECTIVE ANONYMISATION IS NOT GUARANTEED !!!
% The anonymisation implemented here depends on the structure and content
% of the DICOM header and might not be effective in many cases.

if nargin<2
    opts.anonym = 'basic';
end
if ~iscell(hdr)
    hdrout{1} = hdr;
else
    hdrout = hdr;
end

if ~strcmp(opts.anonym,'none')
    for i=1:length(hdrout)
        if isfield(hdrout{i},'PatientName')
            hdrout{i}.PatientName = 'Anonymous';
        end
        if isfield(hdrout{i},'PatientBirthDate')
            t1 = hdrout{i}.PatientBirthDate;
            t2 = hdrout{i}.StudyDate;
            hdrout{i}.PatientAge = round((t2-t1)*10/365.25)/10;
            hdrout{i} = rmfield(hdrout{i},'PatientBirthDate');
        end
        if isfield(hdrout{i},'Filename')
            hdrout{i} = rmfield(hdrout{i},'Filename');
        end
        
        if strcmp(opts.anonym, 'full')
            try, hdrout{i} = rmfield(hdrout{i},'PatientID'); end
            try, hdrout{i} = rmfield(hdrout{i},'PatientSex'); end
            try, hdrout{i} = rmfield(hdrout{i},'PatientAge'); end
            try, hdrout{i} = rmfield(hdrout{i},'PatientSize'); end
            try, hdrout{i} = rmfield(hdrout{i},'PatientWeight'); end
            try, hdrout{i}.CSASeriesHeaderInfo = rmfield(hdrout{i}.CSASeriesHeaderInfo,'UsedPatientWeight'); end
        end
    end
end


%==========================================================================
function tdyhdr = tidy_CSA(csahdr)
% FORMAT tdyhdr = tidy_CSA(csahdr)
% DESCRIPTION:
% To rearrange CSAImageHeaderInfo, CSASeriesHeaderInfo, ... structures
% and make it a simplified, easier to browse structure (Siemens specific).
% USAGE:
% tdyhdr = tidy_CSA(csahdr)
% where csahdr is a structure, content of the CSA field.

tdyhdr = [];
for ccsa = 1:length(csahdr)
    val = get_numaris4_val(csahdr,csahdr(ccsa).name);
    % if val is empty, let's not waste disk space with empty fields
    if ~isempty(val)
        % for elements having a value representation corresponding to a
        % numerical value (or an array of numerical values), convert char
        % into numbers. Here is a list of VR corresponding to numerical
        % values:
        % - DS (Decimal String)
        % - FL (Floating Point Single)
        % - FD (Floating Point Double)
        % - IS (Integer String)
        % - OD (Other Double String)
        % - OF (Other Float String)
        % - SL (Signed Long)
        % - SS (Signed Short)
        % - UL (Unsigned Long)
        % - US (Unsigned Short)
        switch deblank(csahdr(ccsa).vr)
            case {'DS','FL','FD','IS','OD','OF','SL','SS','UL','US'}
                try
                    tmp = zeros(size(val,1),1);
                    for k = 1:size(val,1)
                        tmp(k) = str2num(val(k,:));
                    end
                    val = tmp;
                catch
                    fprintf('Trouble reading CSA header %s (%s) = %s\n',csahdr(ccsa).name, deblank(csahdr(ccsa).vr), val(1,:));
                    val = [];
                end
            otherwise
                val = deblank(val);
        end
        % make sure no "-" in field name (invalid otherwise)
        csahdr(ccsa).name = strrep(csahdr(ccsa).name,'-','');
        tdyhdr.(csahdr(ccsa).name) = val;
    end
end


%==========================================================================
function val = get_numaris4_val(str,name)
% copied from spm_dicom_convert
name = deblank(name);
val  = {};
for i=1:length(str)
    if strcmp(deblank(str(i).name),name)
        if isfield(str(i),'nitems')
            for j=1:str(i).nitems
                if  str(i).item(j).xx(1)
                    val = [val {str(i).item(j).val}];
                end
            end
            %         else
            %             fprintf(1,'Found empty CSA header with name %s\n',name);
        end
        break;
    end
end
val = strvcat(val{:});


%==========================================================================
function ascout = read_ASCII(ascin)
% FORMAT ascout = read_ASCII(ascin)
%
% DESCRIPTION:
% To parse the content of the ASCII part of the DICOM header (Siemens
% specific) and convert it into a matlab structure. Compatible with VA, VB,
% VD and also reading *.SR files (PhoenixProtocols).
%
% USAGE:
% ascout = read_ASCII(ascin)
% where ascin is the char content of either
% hdr.CSASeriesHeaderInfo.MrPhoenixProtocol or
% hdr.CSASeriesHeaderInfo.MrProtocol according to Siemens software version,
% and asc is a matlab structure containing the information contained in the
% ASCCONV BEGIN - ASCCONV END part of the DICOM header.
%
% NOTE: slightly different formatting from read_ascconv in spm_dicom_convert

% only for debugging purpose
ENABLE_DEBUG = false;

% extract the portion between ### ASCCONV BEGIN ### and ### ASCCONV END ###
ascin = regexprep(ascin,'^.*### ASCCONV BEGIN [^#]*###(.*)### ASCCONV END ###.*$','$1');

% split ascin into lines
% asclines = strsplit(ascin,'\n'); % ok for Matlab 2013 and later versions
% replaced by textscan (with delimiter '\n') for compatibility with
% Matlab 2012 and earlier :/...
tmp = textscan(ascin,'%s','delimiter','\n');
asclines = tmp{1};

% do a first cleaning pass over the data:
% - replace indexes [] by () + increment index +1 (C>Matlab indexing)
% - replace double "" and single " by single '
% - replace hexadecimal values (e.g. "0x01") by decimal value
% - delete lines where a field starts or ends by "_"
% - delete end of lines after "#"
asclines = regexprep(asclines,{'\[([0-9]*)\]','["]+','^([^"]*)0x([0-9a-fA-F]*)',' #.*','^.*\._.*$'},{'($1+1)','''','$1hex2dec(''$2'')','',''});

% initialise variables
ascout = [];
clinenum = 1;

while clinenum<length(asclines)
    
    if ENABLE_DEBUG, fprintf('Line #%d - [%s]\n', clinenum, asclines{clinenum}); end
    if ~isempty(asclines{clinenum})
        try
            [tlhrh,tmp] = regexp(asclines{clinenum}, '(?:=)+', 'split', 'match'); %#ok<*ASGLU>
            % first process the name of the parameter (hdrnam)
            % split into fields
            [hdrnam,tmp] = regexp(tlhrh{1}, '(?:\.)+', 'split', 'match');
            % check whether any field starts with a number and if so,
            % replace it by "index<number>..."
            isfirstdigit = cellfun(@(x) isstrprop(x(1),'digit'),hdrnam);
            if any(isfirstdigit)
                for cfield=1:length(hdrnam)
                    if isfirstdigit(cfield)
                        hdrnam{cfield} = ['index' hdrnam{cfield}];
                    end
                end
            end
            % concatenate back the fields
            hdrnam = sprintf('.%s', hdrnam{:});
            % now process the value (hdrval)
            hdrval = strtrim(tlhrh{2});
            % and eval the whole line
            eval(sprintf('ascout%s = %s;', hdrnam, hdrval));
            if ENABLE_DEBUG, fprintf('ascout%s = %s;\n', hdrnam, hdrval); end
        catch
            fprintf('AscConv: Error evaluating [ %s; ]\n', asclines{clinenum});
        end
    end
    clinenum = clinenum+1;
end
