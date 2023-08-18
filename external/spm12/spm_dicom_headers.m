function Headers = spm_dicom_headers(DicomFilenames, Essentials)
% Read header information from DICOM files
% FORMAT Headers = spm_dicom_headers(DicomFilenames [,Essentials])
% DicomFilenames - array of filenames
% Essentials     - if true, then only save the essential parts of the header
%
% Headers        - cell array of headers, one element for each file.
%
% Contents of headers are approximately explained in:
% http://medical.nema.org/standard.html
%
% This code may not work for all cases of DICOM data, as DICOM is an
% extremely complicated "standard".
%__________________________________________________________________________
% Copyright (C) 2002-2018 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_headers.m 7755 2019-12-16 13:19:28Z spm $


DicomFilenames = cellstr(DicomFilenames);

if nargin<2, Essentials = false; end
if ~isa(Essentials,'function_handle')
    if Essentials
        Essentials = @spm_dicom_essentials;
    else
        Essentials = @(x) x;
    end
end

DicomDictionary = load(fullfile(spm('Dir'),'spm_dicom_dict.mat'));

Headers  = {};
if numel(DicomFilenames)>1
    spm_progress_bar('Init',numel(DicomFilenames), ...
        'Reading DICOM headers', 'Files complete');
end
for i=1:numel(DicomFilenames)
    Header = spm_dicom_header(DicomFilenames{i}, DicomDictionary);
    Header = Essentials(Header);
    if ~isempty(Header)
        Headers{end+1} = Header;
    end
    if numel(DicomFilenames)>1, spm_progress_bar('Set',i); end
end
if numel(DicomFilenames)>1, spm_progress_bar('Clear'); end
