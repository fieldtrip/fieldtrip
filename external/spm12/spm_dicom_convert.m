function out = spm_dicom_convert(Headers,opts,RootDirectory,format,OutputDirectory,meta)
% Convert DICOM images into something that SPM can use (e.g. NIfTI)
% FORMAT out = spm_dicom_convert(Headers,opts,RootDirectory,format,OutputDirectory)
% Inputs:
% Headers      - a cell array of DICOM headers from spm_dicom_headers
% opts     - options:
%              'all'      - all DICOM files [default]
%              'mosaic'   - the mosaic images
%              'standard' - standard DICOM files
%              'spect'    - SIEMENS Spectroscopy DICOMs (some formats only)
%                           This will write out a 5D NIFTI containing real
%                           and imaginary part of the spectroscopy time 
%                           points at the position of spectroscopy voxel(s)
%              'raw'      - convert raw FIDs (not implemented)
% RootDirectory - 'flat'       - do not produce file tree [default]
%              With all other options, files will be sorted into
%              directories according to their sequence/protocol names:
%            'date_time'  - Place files under ./<StudyDate-StudyTime>
%            'patid'      - Place files under ./<PatID>
%            'patid_date' - Place files under ./<PatID-StudyDate>
%            'series'     - Place files in series folders, without
%                           creating patient folders
% format   - output format:
%              'nii'      - Single file NIfTI format [default]
%              'img'      - Two file (Headers+img) NIfTI format
%            All images will contain a single 3D dataset, 4D images will
%            not be created.
% OutputDirectory  - output directory name [default: pwd]
% meta     - save metadata as sidecar JSON file [default: false]
%
% Output:
% out      - a struct with a single field .files. out.files contains a
%            cellstring with filenames of created files. If no files are
%            created, a cell with an empty string {''} is returned.
%__________________________________________________________________________
% Copyright (C) 2002-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_convert.m 7714 2019-11-26 11:25:50Z spm $


%-Input parameters
%--------------------------------------------------------------------------
if nargin<2, opts   = 'all';  end
if nargin<3, RootDirectory = 'flat'; end
if nargin<4, format = spm_get_defaults('images.format'); end
if nargin<5, OutputDirectory  = pwd;    end
if nargin<6, meta   = false;  end

%-Select files
%--------------------------------------------------------------------------
[images, other]     = SelectTomographicImages(Headers);
[multiframe,images] = SelectMultiframe(images);
[spect, guff]       = SelectSpectroscopyImages(other);
[mosaic, standard]  = SelectMosaicImages(images);
[standard, guff]    = SelectLastGuff(standard, guff);

if ~isempty(guff)
    warning('spm:dicom','%d files could not be converted from DICOM.', numel(guff));
end

%-Convert files
%--------------------------------------------------------------------------
fmos = {};
fstd = {};
fspe = {};
fmul = {};
if (strcmp(opts,'all') || strcmp(opts,'mosaic')) && ~isempty(mosaic)
    fmos = ConvertMosaic(mosaic,RootDirectory,format,OutputDirectory,meta);
end
if (strcmp(opts,'all') || strcmp(opts,'standard')) && ~isempty(standard)
    fstd = ConvertStandard(standard,RootDirectory,format,OutputDirectory,meta);
end
if (strcmp(opts,'all') || strcmp(opts,'spect')) && ~isempty(spect)
    fspe = ConvertSpectroscopy(spect,RootDirectory,format,OutputDirectory,meta);
end
if (strcmp(opts,'all') || strcmp(opts,'multiframe')) && ~isempty(multiframe)
    fmul = ConvertMultiframes(multiframe,RootDirectory,format,OutputDirectory,meta);
end
out.files = [fmos(:); fstd(:); fspe(:); fmul(:)];
if isempty(out.files)
    out.files = {''};
end


%==========================================================================
% function fnames = ConvertMosaic(Headers,RootDirectory,format,OutputDirectory,meta)
%==========================================================================
function fnames = ConvertMosaic(Headers,RootDirectory,format,OutputDirectory,meta)
spm_progress_bar('Init',length(Headers),'Writing Mosaic', 'Files written');

fnames = cell(length(Headers),1);
for i=1:length(Headers)

    % Output filename
    %----------------------------------------------------------------------
    fnames{i} = getfilelocation(Headers{i}, RootDirectory, 'f', format, OutputDirectory);

    % Image dimensions and data
    %----------------------------------------------------------------------
    nc = Headers{i}.Columns;
    nr = Headers{i}.Rows;

    dim      = [0 0 0];
    dim(3)   = ReadNumberOfImagesInMosaic(Headers{i});
    np       = [nc nr]/ceil(sqrt(dim(3)));
    dim(1:2) = np;
    if ~all(np==floor(np))
        warning('spm:dicom','%s: dimension problem [Num Images=%d, Num Cols=%d, Num Rows=%d].',...
            Headers{i}.Filename,dim(3), nc,nr);
        continue
    end

    % Apparently, this is not the right way of doing it.
    %np = ReadAcquisitionMatrixText(Headers{i});
    %if rem(nc, np(1)) || rem(nr, np(2)),
    %   warning('spm:dicom','%s: %dx%d wont fit into %dx%d.',Headers{i}.Filename,...
    %       np(1), np(2), nc,nr);
    %   return;
    %end;
    %dim    = [np ReadNumberOfImagesInMosaic(Headers{i})];

    mosaic = ReadImageData(Headers{i});
    volume = zeros(dim);
    snnz   = zeros(dim(3), 1);
    for j=1:dim(3)
        img = mosaic((1:np(1))+np(1)*rem(j-1,nc/np(1)), (np(2):-1:1)+np(2)*floor((j-1)/(nc/np(1))));
        snnz(j) = nnz(img) > 0;
        volume(:,:,j) = img;
    end
    d3 = find(snnz, 1, 'last');
    if ~isempty(d3)
        dim(3) = d3;
        volume = volume(:,:,1:dim(3));
    end
    dt  = DetermineDatatype(Headers{1});

    % Orientation information
    %----------------------------------------------------------------------
    % Axial Analyze voxel co-ordinate system:
    % x increases     right to left
    % y increases posterior to anterior
    % z increases  inferior to superior

    % DICOM patient co-ordinate system:
    % x increases     right to left
    % y increases  anterior to posterior
    % z increases  inferior to superior

    % T&T co-ordinate system:
    % x increases      left to right
    % y increases posterior to anterior
    % z increases  inferior to superior

    AnalyzeToDicom = [diag([1 -1 1]) [0 (dim(2)-1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];

    vox    = [Headers{i}.PixelSpacing(:); Headers{i}.SpacingBetweenSlices];
    pos    = Headers{i}.ImagePositionPatient(:);
    orient = reshape(Headers{i}.ImageOrientationPatient,[3 2]);
    orient(:,3) = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end

    % The image position vector is not correct. In dicom this vector points to
    % the upper left corner of the image. Perhaps it is unlucky that this is
    % calculated in the syngo software from the vector pointing to the center of
    % the slice (keep in mind: upper left slice) with the enlarged FoV.
    DicomToPatient = [orient*diag(vox) pos ; 0 0 0 1];
    truepos        = DicomToPatient *[(size(mosaic)-dim(1:2))/2 0 1]';
    DicomToPatient = [orient*diag(vox) truepos(1:3) ; 0 0 0 1];
    PatientToTal   = diag([-1 -1 1 1]);
    mat            = PatientToTal*DicomToPatient*AnalyzeToDicom;


    % Maybe flip the image depending on SliceNormalVector from 0029,1010
    %----------------------------------------------------------------------
    SliceNormalVector = ReadSliceNormalVector(Headers{i});
    if det([reshape(Headers{i}.ImageOrientationPatient,[3 2]) SliceNormalVector(:)])<0
        volume = volume(:,:,end:-1:1);
        mat    = mat*[eye(3) [0 0 -(dim(3)-1)]'; 0 0 0 1];
    end


    % Possibly useful information
    %----------------------------------------------------------------------
    if CheckFields(Headers{i},'AcquisitionTime','MagneticFieldStrength',...
            'MRAcquisitionType','ScanningSequence','RepetitionTime',...
            'EchoTime','FlipAngle','AcquisitionDate')
        tim = datevec(Headers{i}.AcquisitionTime/(24*60*60));
        descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g Mosaic',...
            Headers{i}.MagneticFieldStrength, Headers{i}.MRAcquisitionType,...
            deblank(Headers{i}.ScanningSequence),...
            Headers{i}.RepetitionTime,Headers{i}.EchoTime,Headers{i}.FlipAngle,...
            datestr(Headers{i}.AcquisitionDate),tim(4),tim(5),tim(6));
    else
        descrip = Headers{1}.Modality;
    end

    if ~true % LEFT-HANDED STORAGE
        mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
        volume = flipud(volume);
    end

    % Note that data are no longer scaled by the maximum amount.
    % This may lead to rounding errors in smoothed data, but it
    % will get around other problems.
    RescaleSlope     = 1;
    RescaleIntercept = 0;
    if isfield(Headers{i},'RescaleSlope')
        RescaleSlope = Headers{i}.RescaleSlope;
    end
    if isfield(Headers{i},'RescaleIntercept')
        RescaleIntercept = Headers{i}.RescaleIntercept;
    end
    Nii      = nifti;
    Nii.dat  = file_array(fnames{i},dim,dt,0,RescaleSlope,RescaleIntercept);
    Nii.mat  = mat;
    Nii.mat0 = mat;
    Nii.mat_intent  = 'Scanner';
    Nii.mat0_intent = 'Scanner';
    Nii.descrip     = descrip;
    create(Nii);

    if meta
        Nii = spm_dicom_metadata(Nii,Headers{i});
    end
    
    % Write the data unscaled
    dat           = Nii.dat;
    dat.scl_slope = [];
    dat.scl_inter = [];
    % write out volume at once - see spm_write_plane.m for performance comments
    dat(:,:,:) = volume;
    
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


%==========================================================================
% function fnames = ConvertStandard(Headers,RootDirectory,format,OutputDirectory,meta)
%==========================================================================
function fnames = ConvertStandard(Headers,RootDirectory,format,OutputDirectory,meta)
[Headers,dist] = SortIntoVolumes(Headers);
fnames         = cell(length(Headers),1);
for i=1:length(Headers)
    fnames{i} = WriteVolume(Headers{i},RootDirectory,format,OutputDirectory,dist{i},meta);
end


%==========================================================================
% function [SortedHeaders,dist] = SortIntoVolumes(Headers)
%==========================================================================
function [SortedHeaders,dist] = SortIntoVolumes(Headers)
%
% First of all, sort into volumes based on relevant
% fields in the header.
%

SortedHeaders{1}{1} = Headers{1};
for i=2:length(Headers)
   %orient = reshape(Headers{i}.ImageOrientationPatient,[3 2]);
   %xy1    = Headers{i}.ImagePositionPatient(:)*orient;
    match  = 0;
    if isfield(Headers{i},'CSAImageHeaderInfo') && isfield(Headers{i}.CSAImageHeaderInfo,'name')
        ice1 = sscanf( ...
            strrep(GetNumaris4Val(Headers{i}.CSAImageHeaderInfo,'ICE_Dims'), ...
            'X', '-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
        dimsel = logical([1 1 1 1 1 1 0 0 1]);
    else
        ice1 = [];
    end
    for j=1:length(SortedHeaders)
       %orient = reshape(SortedHeaders{j}{1}.ImageOrientationPatient,[3 2]);
       %xy2    = SortedHeaders{j}{1}.ImagePositionPatient(:)*orient;
        
        % This line is a fudge because of some problematic data that Bogdan,
        % Cynthia and Stefan were trying to convert.  I hope it won't cause
        % problems for others -JA
        % dist2  = sum((xy1-xy2).^2);
        dist2 = 0;
        
        if strcmp(Headers{i}.Modality,'CT') && ...
                strcmp(SortedHeaders{j}{1}.Modality,'CT') % Our CT seems to have shears in slice positions
            dist2 = 0;
        end
        if ~isempty(ice1) && isfield(SortedHeaders{j}{1},'CSAImageHeaderInfo') && numel(SortedHeaders{j}{1}.CSAImageHeaderInfo)>=1 && isfield(SortedHeaders{j}{1}.CSAImageHeaderInfo(1),'name')
            % Replace 'X' in ICE_Dims by '-1'
            ice2 = sscanf( ...
                strrep(GetNumaris4Val(SortedHeaders{j}{1}.CSAImageHeaderInfo,'ICE_Dims'), ...
                'X', '-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
            if ~isempty(ice2)
                identical_ice_dims=all(ice1(dimsel)==ice2(dimsel));
            else
                identical_ice_dims = 0; % have ice1 but not ice2, ->
                % something must be different
            end
        else
            identical_ice_dims = 1; % No way of knowing if there is no CSAImageHeaderInfo
        end
        try
            match = Headers{i}.SeriesNumber            == SortedHeaders{j}{1}.SeriesNumber &&...
                Headers{i}.Rows                        == SortedHeaders{j}{1}.Rows &&...
                Headers{i}.Columns                     == SortedHeaders{j}{1}.Columns &&...
                sum((Headers{i}.ImageOrientationPatient - SortedHeaders{j}{1}.ImageOrientationPatient).^2)<1e-4 &&...
                sum((Headers{i}.PixelSpacing            - SortedHeaders{j}{1}.PixelSpacing).^2)<1e-4 && ...
                identical_ice_dims && dist2<1e-3;
            %if (Headers{i}.AcquisitionNumber ~= Headers{i}.InstanceNumber) || ...
            %   (SortedHeaders{j}{1}.AcquisitionNumber ~= SortedHeaders{j}{1}.InstanceNumber)
            %    match = match && (Headers{i}.AcquisitionNumber == SortedHeaders{j}{1}.AcquisitionNumber)
            %end
            % For raw image data, tell apart real/complex or phase/magnitude
            if isfield(Headers{i},'ImageType') && isfield(SortedHeaders{j}{1}, 'ImageType')
                match = match && strcmp(Headers{i}.ImageType, SortedHeaders{j}{1}.ImageType);
            end
            if isfield(Headers{i},'SequenceName') && isfield(SortedHeaders{j}{1}, 'SequenceName')
                match = match && strcmp(Headers{i}.SequenceName, SortedHeaders{j}{1}.SequenceName);
            end
            if isfield(Headers{i},'SeriesInstanceUID') && isfield(SortedHeaders{j}{1}, 'SeriesInstanceUID')
                match = match && strcmp(Headers{i}.SeriesInstanceUID, SortedHeaders{j}{1}.SeriesInstanceUID);
            end
            if isfield(Headers{i},'EchoNumbers')  && isfield(SortedHeaders{j}{1}, 'EchoNumbers')
                match = match && Headers{i}.EchoNumbers == SortedHeaders{j}{1}.EchoNumbers;
            end
            if isfield(Headers{i},    'GE_ImageType') && numel(   Headers{i}.GE_ImageType)==1 && ...
               isfield(SortedHeaders{j}{1}, 'GE_ImageType') && numel(SortedHeaders{j}{1}.GE_ImageType)==1
                match = match && Headers{i}.GE_ImageType == SortedHeaders{j}{1}.GE_ImageType;
            end
        catch
            match = 0;
        end
        if match
            SortedHeaders{j}{end+1} = Headers{i};
            break;
        end
    end
    if ~match
        SortedHeaders{end+1}{1} = Headers{i};
    end
end

%
% Secondly, sort volumes into ascending/descending
% slices depending on .ImageOrientationPatient field.
%

SortedHeaders2 = {};
for j=1:length(SortedHeaders)
    orient = reshape(SortedHeaders{j}{1}.ImageOrientationPatient,[3 2]);
    proj   = null(orient');
    if det([orient proj])<0, proj = -proj; end

    z      = zeros(length(SortedHeaders{j}),1);
    for i=1:length(SortedHeaders{j})
        z(i)  = SortedHeaders{j}{i}.ImagePositionPatient(:)'*proj;
    end
    [z,index]        = sort(z);
    SortedHeaders{j} = SortedHeaders{j}(index);
    if length(SortedHeaders{j})>1
        % dist      = diff(z);
        if any(diff(z)==0)
            tmp = SortIntoVolumesAgain(SortedHeaders{j});
            SortedHeaders{j} = tmp{1};
            SortedHeaders2 = {SortedHeaders2{:} tmp{2:end}};
        end
    end
end
SortedHeaders = {SortedHeaders{:} SortedHeaders2{:}};
dist          = cell(length(SortedHeaders),1);
for j=1:length(SortedHeaders)
    if length(SortedHeaders{j})>1
        orient = reshape(SortedHeaders{j}{1}.ImageOrientationPatient,[3 2]);
        proj   = null(orient');
        if det([orient proj])<0, proj = -proj; end
        z      = zeros(length(SortedHeaders{j}),1);
        for i=1:length(SortedHeaders{j})
            z(i)  = SortedHeaders{j}{i}.ImagePositionPatient(:)'*proj;
        end
        dist{j} = diff(sort(z));
        if sum((dist{j}-mean(dist{j})).^2)/length(dist{j})>1e-4
            fprintf('***************************************************\n');
            fprintf('* VARIABLE SLICE SPACING                          *\n');
            fprintf('* This may be due to missing DICOM files.         *\n');
            PatientID = 'anon';
            if CheckFields(SortedHeaders{j}{1}, 'PatientID'), PatientID = deblank(SortedHeaders{j}{1}.PatientID); end
            if CheckFields(SortedHeaders{j}{1}, 'SeriesNumber', 'AcquisitionNumber', 'InstanceNumber')
                fprintf('*    %s / %d / %d / %d \n',...
                    PatientID, SortedHeaders{j}{1}.SeriesNumber, ...
                    SortedHeaders{j}{1}.AcquisitionNumber, SortedHeaders{j}{1}.InstanceNumber);
                fprintf('*                                                 *\n');
            end
            fprintf('*  %20.4g                           *\n', dist{j});
            fprintf('***************************************************\n');
        end
    end
end


%==========================================================================
% function SortedHeaders = SortIntoVolumesAgain(Headers)
%==========================================================================
function SortedHeaders = SortIntoVolumesAgain(Headers)
if ~isfield(Headers{1},'InstanceNumber')
    fprintf('***************************************************\n');
    fprintf('* The slices may be all mixed up and the data     *\n');
    fprintf('* not really usable.  Talk to your physicists     *\n');
    fprintf('* about this.                                     *\n');
    fprintf('***************************************************\n');
    SortedHeaders = {Headers};
    return;
end

fprintf('***************************************************\n');
fprintf('* The AcquisitionNumber counter does not appear   *\n');
fprintf('* to be changing from one volume to another.      *\n');
fprintf('* Another possible explanation is that the same   *\n');
fprintf('* DICOM slices are used multiple times.           *\n');
%fprintf('* Talk to your MR sequence developers or scanner  *\n');
%fprintf('* supplier to have this fixed.                    *\n');
fprintf('* The conversion is having to guess how slices    *\n');
fprintf('* should be arranged into volumes.                *\n');
PatientID = 'anon';
if CheckFields(Headers{1},'PatientID'), PatientID = deblank(Headers{1}.PatientID); end
if CheckFields(Headers{1},'SeriesNumber','AcquisitionNumber')
    fprintf('*    %s / %d / %d\n',...
        PatientID, Headers{1}.SeriesNumber, ...
        Headers{1}.AcquisitionNumber);
end
fprintf('***************************************************\n');

z      = zeros(length(Headers),1);
t      = zeros(length(Headers),1);
d      = zeros(length(Headers),1);
orient = reshape(Headers{1}.ImageOrientationPatient,[3 2]);
proj   = null(orient');
if det([orient proj])<0, proj = -proj; end

for i=1:length(Headers)
    z(i)  = Headers{i}.ImagePositionPatient(:)'*proj;
    t(i)  = Headers{i}.InstanceNumber;
end
% msg = 0;
[t,index]      = sort(t);
Headers = Headers(index);
z              = z(index);
msk            = find(diff(t)==0);
if any(msk)
    % fprintf('***************************************************\n');
    % fprintf('* These files have the same InstanceNumber:       *\n');
    % for i=1:length(msk),
    %    [tmp,nam1,ext1] = fileparts(Headers{msk(i)}.Filename);
    %    [tmp,nam2,ext2] = fileparts(Headers{msk(i)+1}.Filename);
    %    fprintf('* %s%s = %s%s (%d)\n', nam1,ext1,nam2,ext2, Headers{msk(i)}.InstanceNumber);
    % end;
    % fprintf('***************************************************\n');
    index = [true ; diff(t)~=0];
    t     = t(index);
    z     = z(index);
    d     = d(index);
    Headers  = Headers(index);
end

%if any(diff(sort(t))~=1), msg = 1; end;
[z,index]     = sort(z);
Headers       = Headers(index);
t             = t(index);
SortedHeaders = {};
while ~all(d)
    i  = find(~d);
    i  = i(1);
    i  = find(z==z(i));
    [t(i),si] = sort(t(i));
    Headers(i)   = Headers(i(si));
    for i1=1:length(i)
        if length(SortedHeaders)<i1, SortedHeaders{i1} = {}; end
        SortedHeaders{i1} = {SortedHeaders{i1}{:} Headers{i(i1)}};
    end
    d(i) = 1;
end

msg = 0;
len = length(SortedHeaders{1});
for i=2:length(SortedHeaders)
    if length(SortedHeaders{i}) ~= len
        msg = 1;
        break;
    end
end
if msg
    fprintf('***************************************************\n');
    fprintf('* There are missing DICOM files, so the the       *\n');
    fprintf('* resulting volumes may be messed up.             *\n');
    PatientID = 'anon';
    if CheckFields(Headers{1},'PatientID'), PatientID = deblank(Headers{1}.PatientID); end
    if CheckFields(Headers{1},'SeriesNumber','AcquisitionNumber')
        fprintf('*    %s / %d / %d\n',...
            PatientID, Headers{1}.SeriesNumber, ...
            Headers{1}.AcquisitionNumber);
    end
    fprintf('***************************************************\n');
end


%==========================================================================
% function fname = WriteVolume(Headers, RootDirectory, format, OutputDirectory, dist, meta)
%==========================================================================
function fname = WriteVolume(Headers, RootDirectory, format, OutputDirectory, dist, meta)

% Output filename
%--------------------------------------------------------------------------
fname = getfilelocation(Headers{1}, RootDirectory,'s',format,OutputDirectory);

% Image dimensions
%--------------------------------------------------------------------------
nc = Headers{1}.Columns;
nr = Headers{1}.Rows;

if length(Headers) == 1 && isfield(Headers{1},'NumberOfFrames') && Headers{1}.NumberOfFrames > 1
    if isfield(Headers{1},'ImagePositionPatient') &&...
       isfield(Headers{1},'ImageOrientationPatient') &&...
       isfield(Headers{1},'SliceThickness') &&...
       isfield(Headers{1},'StartOfPixelData') &&...
       isfield(Headers{1},'SizeOfPixelData')

       orient           = reshape(Headers{1}.ImageOrientationPatient,[3 2]);
       orient(:,3)      = null(orient');
       if det(orient)<0, orient(:,3) = -orient(:,3); end
       slicevec         = orient(:,3);
   
       Headers_temp = cell(1,Headers{1}.NumberOfFrames); % alternative: NumberofSlices
       Headers_temp{1} = Headers{1};
       Headers_temp{1}.SizeOfPixelData           = Headers{1}.SizeOfPixelData / Headers{1}.NumberOfFrames;
       for sn = 2 : Headers{1}.NumberOfFrames
           Headers_temp{sn}                      = Headers{1};
           Headers_temp{sn}.ImagePositionPatient = Headers{1}.ImagePositionPatient + (sn-1) * Headers{1}.SliceThickness * slicevec;
           Headers_temp{sn}.SizeOfPixelData      = Headers_temp{1}.SizeOfPixelData;
           Headers_temp{sn}.StartOfPixelData     = Headers{1}.StartOfPixelData + (sn-1) * Headers_temp{1}.SizeOfPixelData;
       end
       Headers = Headers_temp;
   else
       error('spm_dicom_convert:WriteVolume','TAGS missing in DICOM file.');
   end
end

dim    = [nc nr length(Headers)];
dt     = DetermineDatatype(Headers{1});

% Orientation information
%--------------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

AnalyzeToDicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
PatientToTal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions

R  = [reshape(Headers{1}.ImageOrientationPatient,3,2)*diag(Headers{1}.PixelSpacing); 0 0];
x1 = [1;1;1;1];
y1 = [Headers{1}.ImagePositionPatient(:); 1];

if length(Headers)>1
    x2 = [1;1;dim(3); 1];
    y2 = [Headers{end}.ImagePositionPatient(:); 1];
else
    orient           = reshape(Headers{1}.ImageOrientationPatient,[3 2]);
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
    if CheckFields(Headers{1},'SliceThickness')
        z = Headers{1}.SliceThickness;
    else
        z = 1;
    end
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
DicomToPatient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat            = PatientToTal*DicomToPatient*AnalyzeToDicom;

% Possibly useful information
%--------------------------------------------------------------------------
if CheckFields(Headers{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate')
    if isfield(Headers{1},'ScanOptions')
        ScanOptions = Headers{1}.ScanOptions;
    else
        ScanOptions = 'no';
    end
    tim = datevec(Headers{1}.AcquisitionTime/(24*60*60));
    descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg/SO=%s %s %d:%d:%.5g',...
        Headers{1}.MagneticFieldStrength, Headers{1}.MRAcquisitionType,...
        deblank(Headers{1}.ScanningSequence),...
        Headers{1}.RepetitionTime,Headers{1}.EchoTime,Headers{1}.FlipAngle,...
        ScanOptions,...
        datestr(Headers{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
    descrip = Headers{1}.Modality;
end

if ~true % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end

% Write the image volume
%--------------------------------------------------------------------------
spm_progress_bar('Init',length(Headers),['Writing ' fname], 'Planes written');
pinfos = [ones(length(Headers),1) zeros(length(Headers),1)];
for i=1:length(Headers)
    if isfield(Headers{i},'RescaleSlope'),      pinfos(i,1) = Headers{i}.RescaleSlope;      end 
    if isfield(Headers{i},'RescaleIntercept'),  pinfos(i,2) = Headers{i}.RescaleIntercept;  end

    % Philips do things differently. The following is for using their scales instead.
    %     Chenevert, Thomas L., et al. "Errors in quantitative image analysis due to
    %     platform-dependent image scaling." Translational oncology 7.1 (2014): 65-71.
    if isfield(Headers{i},'MRScaleSlope'), pinfos(i,1)     = 1/Headers{i}.MRScaleSlope;                 end
    if isfield(Headers{i},'MRScaleIntercept'), pinfos(i,2) =  -Headers{i}.MRScaleIntercept*pinfos(i,1); end

end

if any(any(diff(pinfos,1)))
    % Ensure random numbers are reproducible (see later)
    % when intensities are dithered to prevent aliasing effects.
    rand('state',0);
end

volume = zeros(dim);
for i=1:length(Headers)
    plane = ReadImageData(Headers{i});

    if any(any(diff(pinfos,1)))
        % This is to prevent aliasing effects in any subsequent histograms
        % of the data (eg for mutual information coregistration).
        % It's a bit inelegant, but probably necessary for when slices are
        % individually rescaled.
        plane = double(plane) + rand(size(plane)) - 0.5;
    end

    if pinfos(i,1)~=1, plane = plane*pinfos(i,1); end
    if pinfos(i,2)~=0, plane = plane+pinfos(i,2); end

    plane = fliplr(plane);
    if ~true, plane = flipud(plane); end % LEFT-HANDED STORAGE
    volume(:,:,i) = plane;
    spm_progress_bar('Set',i);
end

if ~any(any(diff(pinfos,1)))
    % Same slopes and intercepts for all slices
    pinfo = pinfos(1,:);
else
    % Variable slopes and intercept (maybe PET/SPECT)
    mx = max(volume(:));
    mn = min(volume(:));

    %  Slope and Intercept
    %  32767*pinfo(1) + pinfo(2) = mx
    % -32768*pinfo(1) + pinfo(2) = mn
    % pinfo = ([32767 1; -32768 1]\[mx; mn])';

    % Slope only
    dt    = 'int16-be';
    pinfo = [max(mx/32767,-mn/32768) 0];
end

Nii      = nifti;
Nii.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
Nii.mat  = mat;
Nii.mat0 = mat;
Nii.mat_intent  = 'Scanner';
Nii.mat0_intent = 'Scanner';
Nii.descrip     = descrip;
create(Nii);

if meta
    Nii = spm_dicom_metadata(Nii, Headers{1});
end

Nii.dat(:,:,:) = volume;
spm_progress_bar('Clear');

if sum((dist-mean(dist)).^2)/length(dist)>1e-4
    % Adjusting for variable slice thickness
    % This is sometimes required for CT data because these scans are often 
    % acquired with a gantry tilt and thinner slices near the brain stem. 
    % If uncorrected, the resulting NIfTI image can appear distorted. Here, 
    % the image is resampled to compensate for this effect.
    %----------------------------------------------------------------------
    
    fprintf('***************************************************\n');
    fprintf('* Adjusting for variable slice thickness.         *\n');
    fprintf('***************************************************\n');

    ovx = sqrt(sum(mat(1:3,1:3).^2));  
    nvx = [ovx(1:2) max(min(dist),1)]; % Set new slice thickness in thick-slice direction to minimum of dist
    
    Nii = nifti(fname);
    img = Nii.dat(:,:,:);
    d   = size(img);    
    
    csdist = cumsum(dist);  
    csdist = [1; 1 + csdist];
    
    x       = zeros([d(1:3) 3],'single');
    [x1,x2] = meshgrid(single(1:d(2)),single(1:d(1)));
    for i=1:d(3)
        x(:,:,i,1) = x1;
        x(:,:,i,2) = x2;
        x(:,:,i,3) = csdist(i);
    end

    X = x(:,:,:,1);
    Y = x(:,:,:,2);
    Z = x(:,:,:,3);
    clear x
    
    if ~(csdist(end) == floor(csdist(end)))
        [Xq,Yq,Zq] = meshgrid(single(1:d(2)),single(1:d(1)),single([1:nvx(3):csdist(end) csdist(end)]));
    else
        [Xq,Yq,Zq] = meshgrid(single(1:d(2)),single(1:d(1)),single(1:nvx(3):csdist(end)));
    end

    img = interp3(X,Y,Z,img,Xq,Yq,Zq,'linear'); 
    
    ndim = size(img);   
    
    % Adjust orientation matrix
    D   = diag([ovx./nvx 1]);
    mat = mat/D;    
    
    N             = nifti;
    N.dat         = file_array(fname,ndim,dt,0,pinfo(1),pinfo(2));
    N.mat         = mat;
    N.mat0        = mat;
    N.mat_intent  = 'Scanner';
    N.mat0_intent = 'Scanner';
    N.descrip     = descrip;
    create(N);
    N.dat(:,:,:) = img;
end

%==========================================================================
% function fnames = ConvertSpectroscopy(Headers, RootDirectory, format, OutputDirectory, meta)
%==========================================================================
function fnames = ConvertSpectroscopy(Headers, RootDirectory, format, OutputDirectory, meta)
fnames = cell(length(Headers),1);
for i=1:length(Headers)
    fnames{i} = WriteSpectroscopyVolume(Headers(i), RootDirectory, format, OutputDirectory, meta);
end


%==========================================================================
% function fname = WriteSpectroscopyVolume(Headers, RootDirectory, format, OutputDirectory, meta)
%==========================================================================
function fname = WriteSpectroscopyVolume(Headers, RootDirectory, format, OutputDirectory, meta)
% Output filename
%-------------------------------------------------------------------
fname = getfilelocation(Headers{1}, RootDirectory,'S',format,OutputDirectory);

% private field to use - depends on SIEMENS software version
if isfield(Headers{1}, 'CSANonImageHeaderInfoVA')
    privdat = Headers{1}.CSANonImageHeaderInfoVA;
elseif isfield(Headers{1}, 'CSANonImageHeaderInfoVB')
    privdat = Headers{1}.CSANonImageHeaderInfoVB;
else
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end

% Image dimensions
%--------------------------------------------------------------------------
nc = GetNumaris4NumVal(privdat,'Columns');
nr = GetNumaris4NumVal(privdat,'Rows');
% Guess number of timepoints in file - I don't know for sure whether this should be
% 'DataPointRows'-by-'DataPointColumns', 'SpectroscopyAcquisitionDataColumns'
% or sSpecPara.lVectorSize from SIEMENS ASCII header
% ntp = GetNumaris4NumVal(privdat,'DataPointRows')*GetNumaris4NumVal(privdat,'DataPointColumns');
ac = ReadAscconv(Headers{1});
try
    ntp = ac.sSpecPara.lVectorSize;
catch
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end
dim    = [nc nr numel(Headers) 2 ntp];
dt     = spm_type('float32'); % Fixed datatype

% Orientation information
%--------------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

AnalyzeToDicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
PatientToTal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions
shift_vx       = [eye(4,3) [.5; .5; 0; 1]];

orient         = reshape(GetNumaris4NumVal(privdat, 'ImageOrientationPatient'), [3 2]);
ps             = GetNumaris4NumVal(privdat, 'PixelSpacing');
if nc*nr == 1
    % Single Voxel Spectroscopy (based on the following information from SIEMENS)
    %----------------------------------------------------------------------
    % NOTE: Internally the position vector of the CSI matrix shows to the
    % outer border of the first voxel. Therefore the position vector has to
    % be corrected. (Note: The convention of Siemens spectroscopy raw data
    % is in contrast to the DICOM standard where the position vector points
    % to the center of the first voxel.)
    %----------------------------------------------------------------------
    % SIEMENS decides which definition to use based on the contents of the
    % 'PixelSpacing' internal header field. If it has non-zero values,
    % assume DICOM convention. If any value is zero, assume SIEMENS
    % internal convention for this direction.
    % Note that in SIEMENS code, there is a shift when PixelSpacing is
    % zero. Here, the shift seems to be necessary when PixelSpacing is
    % non-zero. This may indicate more fundamental problems with
    % orientation decoding.
    if ps(1) == 0 % row
        ps(1) = GetNumaris4NumVal(privdat, 'VoiPhaseFoV');
        shift_vx(1,4) = 0;
    end
    if ps(2) == 0 % col
        ps(2) = GetNumaris4NumVal(privdat, 'VoiReadoutFoV');
        shift_vx(2,4) = 0;
    end
end
pos = GetNumaris4NumVal(privdat, 'ImagePositionPatient');
% for some reason, pixel spacing needs to be swapped
R  = [orient*diag(ps([2 1])); 0 0];
x1 = [1;1;1;1];
y1 = [pos; 1];

if length(Headers)>1
    error('spm_dicom_convert:spectroscopy',...
        'Don''t know how to handle multislice spectroscopy data.');
else
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
    z = GetNumaris4NumVal(privdat,...
        'VoiThickness');
    if isempty(z)
        z = GetNumaris4NumVal(privdat,...
            'SliceThickness');
    end
    if isempty(z)
        warning('spm_dicom_convert:spectroscopy',...
            'Can not determine voxel thickness.');
        z = 1;
    end
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
DicomToPatient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat              = patient_to_tal*DicomToPatient*shift_vx*AnalyzeToDicom;

% Possibly useful information
%--------------------------------------------------------------------------
if CheckFields(Headers{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate')
    tim = datevec(Headers{1}.AcquisitionTime/(24*60*60));
    descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g',...
        Headers{1}.MagneticFieldStrength, Headers{1}.MRAcquisitionType,...
        deblank(Headers{1}.ScanningSequence),...
        Headers{1}.RepetitionTime,Headers{1}.EchoTime,Headers{1}.FlipAngle,...
        datestr(Headers{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
    descrip = Headers{1}.Modality;
end

if ~true % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end

% Write the image volume
%--------------------------------------------------------------------------
Nii    = nifti;
pinfo  = [1 0];
if isfield(Headers{1},'RescaleSlope'),      pinfo(1) = Headers{1}.RescaleSlope;     end
if isfield(Headers{1},'RescaleIntercept'),  pinfo(2) = Headers{1}.RescaleIntercept; end

% Philips do things differently. The following is for using their scales instead.
%     Chenevert, Thomas L., et al. "Errors in quantitative image analysis due to
%     platform-dependent image scaling." Translational oncology 7.1 (2014): 65-71.
if isfield(Headers{1},'MRScaleSlope'),     pinfo(1) = 1/Headers{1}.MRScaleSlope;              end
if isfield(Headers{1},'MRScaleIntercept'), pinfo(2) =  -Headers{1}.MRScaleIntercept*pinfo(1); end


Nii.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
Nii.mat  = mat;
Nii.mat0 = mat;
Nii.mat_intent  = 'Scanner';
Nii.mat0_intent = 'Scanner';
Nii.descrip     = descrip;
% Store LCMODEL control/raw info in Nii.extras
Nii.extras      = struct('MagneticFieldStrength', GetNumaris4NumVal(privdat,'MagneticFieldStrength'),...
                         'TransmitterReferenceAmplitude', GetNumaris4NumVal(privdat,'TransmitterReferenceAmplitude'),...
                         'ImagingFrequency', GetNumaris4NumVal(privdat,'ImagingFrequency'),...
                         'EchoTime', GetNumaris4NumVal(privdat,'EchoTime'),...
                         'RealDwellTime', GetNumaris4NumVal(privdat,'RealDwellTime'));
create(Nii);

if meta
    Nii = spm_dicom_metadata(Nii,Headers{1});
end

% Read data, swap dimensions
data = permute(reshape(ReadSpectData(Headers{1},ntp),dim([4 5 1 2 3])), ...
                [3 4 5 1 2]);
% plane = fliplr(plane);

Nii.dat(:,:,:,:,:) = data;


%==========================================================================
% function [images,guff] = SelectTomographicImages(Headers)
%==========================================================================
function [images,guff] = SelectTomographicImages(Headers)
images = {};
guff   = {};
for i=1:length(Headers)
    if ~CheckFields(Headers{i},'Modality') || ...
            ~(strcmp(Headers{i}.Modality,'MR') || ...
              strcmp(Headers{i}.Modality,'PT') || ...
              strcmp(Headers{i}.Modality,'NM') || ...
              strcmp(Headers{i}.Modality,'CT'))
        if CheckFields(Headers{i},'Modality')
            fprintf('File "%s" can not be converted because it is of type "%s", which is not MRI, CT, NM or PET.\n', ...
                    Headers{i}.Filename, Headers{i}.Modality);
        else
            fprintf('File "%s" can not be converted because it does not encode an image.\n', Headers{i}.Filename);
        end
        guff = [guff(:)',Headers(i)];
 
    elseif ~CheckFields(Headers{i},'StartOfPixelData','SamplesPerPixel',...
            'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation')
        fprintf('Cant find "Image Pixel" information for "%s".\n',Headers{i}.Filename);
        guff = [guff(:)',Headers(i)];
        
    elseif ~(CheckFields(Headers{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient') ...
            || isfield(Headers{i},'CSANonImageHeaderInfoVA') || isfield(Headers{i},'CSANonImageHeaderInfoVB'))
        if isfield(Headers{i},'SharedFunctionalGroupsSequence') || isfield(Headers{i},'PerFrameFunctionalGroupsSequence')
            fprintf(['\n"%s" appears to be multi-frame DICOM.\n'...
                     'Converting these data is still experimental and has only been tested on a very small number of\n'...
                     'multiframe DICOM files. Feedback about problems would be appreciated - particularly if you can\n'...
                     'give us examples of problematic data (providing there are no subject confidentiality issues).\n\n'], Headers{i}.Filename);
            images = [images(:)',Headers(i)];
        else
            fprintf('No "Image Plane" information for "%s".\n',Headers{i}.Filename);
            guff = [guff(:)',Headers(i)];
        end

    elseif ~CheckFields(Headers{i},'SeriesNumber','AcquisitionNumber','InstanceNumber')
       %disp(['Cant find suitable filename info for "' Headers{i}.Filename '".']);
        if ~isfield(Headers{i},'SeriesNumber')
            fprintf('Setting SeriesNumber to 1.\n');
            Headers{i}.SeriesNumber = 1;
            images = [images(:)',Headers(i)];
        end
        if ~isfield(Headers{i},'AcquisitionNumber')
            if isfield(Headers{i},'Manufacturer') && ~isempty(strfind(upper(Headers{1}.Manufacturer), 'PHILIPS'))
                % Philips oddity
                if isfield(Headers{i},'InstanceNumber')
                    Headers{i}.AcquisitionNumber = Headers{i}.InstanceNumber;
                else
                    fprintf('Setting AcquisitionNumber to 1.\n');
                    Headers{i}.AcquisitionNumber = 1;
                end
            else
                fprintf('Setting AcquisitionNumber to 1.\n');
                Headers{i}.AcquisitionNumber = 1;
            end
            images = [images(:)',Headers(i)];
        end
        if ~isfield(Headers{i},'InstanceNumber')
            fprintf('Setting InstanceNumber to 1.\n');
            Headers{i}.InstanceNumber = 1;
            images = [images(:)',Headers(i)];
        end
    %elseif isfield(Headers{i},'Private_2001_105f'),
    %    % This field corresponds to: > Stack Sequence 2001,105F SQ VNAP, COPY
    %    % http://www.medical.philips.com/main/company/connectivity/mri/index.html
    %    % No documentation about this private field is yet available.
    %    disp('Cant yet convert Philips Intera DICOM.');
    %    guff = {guff{:},Headers{i}};
    else
        images = [images(:)',Headers(i)];
    end
end


%==========================================================================
% function [multiframe,other] = SelectMultiframe(Headers)
%==========================================================================
function [multiframe,other] = SelectMultiframe(Headers)
multiframe = {};
other      = {};
for i=1:length(Headers)
    if isfield(Headers{i},'SharedFunctionalGroupsSequence') || isfield(Headers{i},'PerFrameFunctionalGroupsSequence')
        multiframe = [multiframe(:)',Headers(i)];
    else
        other      = [other(:)',Headers(i)];
    end
end


%==========================================================================
% function [mosaic,standard] = SelectMosaicImages(Headers)
%==========================================================================
function [mosaic,standard] = SelectMosaicImages(Headers)
mosaic   = {};
standard = {};
for i=1:length(Headers)
    if ~CheckFields(Headers{i},'ImageType','CSAImageHeaderInfo') ||...
            isfield(Headers{i}.CSAImageHeaderInfo,'junk') ||...
            isempty(ReadAcquisitionMatrixText(Headers{i})) ||...
            isempty(ReadNumberOfImagesInMosaic(Headers{i})) ||...
            ReadNumberOfImagesInMosaic(Headers{i}) == 0
        % NumberOfImagesInMosaic seems to be set to zero for pseudo images
        % containing e.g. online-fMRI design matrices, don't treat them as
        % mosaics
        standard = [standard, Headers(i)];
    else
        mosaic   = [mosaic,   Headers(i)];
    end
end


%==========================================================================
% function [spect,images] = SelectSpectroscopyImages(Headers)
%==========================================================================
function [spect,images] = SelectSpectroscopyImages(Headers)
spectsel = false(1,numel(Headers));
for i=1:numel(Headers)
    if isfield(Headers{i},'SOPClassUID')
        spectsel(i) = strcmp(Headers{i}.SOPClassUID,'1.3.12.2.1107.5.9.1');
    end
end
spect  = Headers(spectsel);
images = Headers(~spectsel);


%==========================================================================
% function [standard, guff] = SelectLastGuff(standard, guff)
%==========================================================================
function [standard, guff] = SelectLastGuff(standard, guff)
% See https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;5b69d495.1108
i = find(cellfun(@(x) ~isfield(x,'ImageOrientationPatient'),standard));
guff = [guff, standard(i)];
standard(i) = [];


%==========================================================================
% function ok = CheckFields(Headers,varargin)
%==========================================================================
function ok = CheckFields(Headers,varargin)
ok = 1;
for i=1:(nargin-1)
    if ~isfield(Headers,varargin{i})
        ok = 0;
        break;
    end
end


%==========================================================================
% function clean = StripUnwantedChars(dirty)
%==========================================================================
function clean = StripUnwantedChars(dirty)
msk = (dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
      (dirty>='0'&dirty<='9') | dirty=='_';
clean = dirty(msk);


%==========================================================================
% function img = ReadImageData(Header)
%==========================================================================
function img = ReadImageData(Header)
img = [];

if Header.SamplesPerPixel ~= 1
    warning('spm:dicom','%s: SamplesPerPixel = %d - cant be an MRI.', Header.Filename, Header.SamplesPerPixel);
    return;
end

prec = ['ubit' num2str(Header.BitsAllocated) '=>' 'uint32'];

if isfield(Header,'TransferSyntaxUID') && strcmp(Header.TransferSyntaxUID,'1.2.840.10008.1.2.2') && strcmp(Header.VROfPixelData,'OW')
    fp = fopen(Header.Filename,'r','ieee-be');
else
    fp = fopen(Header.Filename,'r','ieee-le');
end
if fp==-1
    warning('spm:dicom','%s: Cant open file.', Header.Filename);
    return;
end

if isfield(Header,'NumberOfFrames')
    NFrames = Header.NumberOfFrames;
else
    NFrames = 1;
end

if isfield(Header,'TransferSyntaxUID')
    switch(Header.TransferSyntaxUID)
    case {'1.2.840.10008.1.2.4.50','1.2.840.10008.1.2.4.51',... % 8 bit JPEG & 12 bit JPEG
          '1.2.840.10008.1.2.4.57','1.2.840.10008.1.2.4.70',... % lossless NH JPEG & lossless NH, 1st order
          '1.2.840.10008.1.2.4.80','1.2.840.10008.1.2.4.81',... % lossless JPEG-LS & near lossless JPEG-LS
          '1.2.840.10008.1.2.4.90','1.2.840.10008.1.2.4.91',... % lossless JPEG 2000 & possibly lossy JPEG 2000, Part 1
          '1.2.840.10008.1.2.4.92','1.2.840.10008.1.2.4.93' ... % lossless JPEG 2000 & possibly lossy JPEG 2000, Part 2
         }
        % try to read PixelData as JPEG image
        fseek(fp,Header.StartOfPixelData,'bof');
        fread(fp,2,'uint16'); % uint16 encoding 65534/57344 (Item)
        offset = double(fread(fp,1,'uint32')); % followed by 4 0 0 0
        fread(fp,2,'uint16'); % uint16 encoding 65534/57344 (Item)
        fread(fp,offset);
        
        sz  = double(fread(fp,1,'*uint32'));
        img = fread(fp,sz,'*uint8');

        % Next uint16 seem to encode 65534/57565 (SequenceDelimitationItem), followed by 0 0

        % save PixelData into temp file - imread and its subroutines can only
        % read from file, not from memory
        tfile = tempname;
        tfp   = fopen(tfile,'w+');
        fwrite(tfp,img,'uint8');
        fclose(tfp);

        % read decompressed data, transpose to match DICOM row/column order
        img = uint32(imread(tfile)');
        delete(tfile);
    case {'1.2.840.10008.1.2.4.94' ,'1.2.840.10008.1.2.4.95' ,... % JPIP References & JPIP Referenced Deflate Transfer
          '1.2.840.10008.1.2.4.100','1.2.840.10008.1.2.4.101',... % MPEG2 MP@ML & MPEG2 MP@HL
          '1.2.840.10008.1.2.4.102',                          ... % MPEG-4 AVC/H.264 High Profile and BD-compatible
         }
         warning('spm:dicom',[Header.Filename ': cant deal with JPIP/MPEG data (' Header.TransferSyntaxUID ')']);
    otherwise
        fseek(fp,Header.StartOfPixelData,'bof');
        img = fread(fp,Header.Rows*Header.Columns*NFrames,prec);
    end
else
    fseek(fp,Header.StartOfPixelData,'bof');
    img = fread(fp,Header.Rows*Header.Columns*NFrames,prec);
end
fclose(fp);
if numel(img)~=Header.Rows*Header.Columns*NFrames
    error([Header.Filename ': cant read whole image']);
end

if Header.PixelRepresentation
    % Signed data - done this way because bitshift only
    % works with signed data.  Negative values are stored
    % as 2s complement.
    neg      = logical(bitshift(bitand(img,uint32(2^Header.HighBit)), -Header.HighBit));
    msk      = (2^Header.HighBit - 1);
    img      = double(bitand(img,msk));
    img(neg) = img(neg)-2^(Header.HighBit);
else
    % Unsigned data
    msk      = (2^(Header.HighBit+1) - 1);
    img      = double(bitand(img,msk));
end

img = reshape(img,[Header.Columns,Header.Rows,NFrames]);


%==========================================================================
% function img = ReadSpectData(Header,privdat)
%==========================================================================
function img = ReadSpectData(Header,ntp)
% Data is stored as complex float32 values, timepoint by timepoint, voxel
% by voxel. Reshaping is done in WriteSpectroscopyVolume.
if ntp*2*4 ~= Header.SizeOfCSAData
    warning('spm:dicom', [Header.Filename,': Data size mismatch.']);
end
fp = fopen(Header.Filename,'r','ieee-le');
fseek(fp,Header.StartOfCSAData,'bof');
img = fread(fp,2*ntp,'float32');
fclose(fp);


%==========================================================================
% function nrm = ReadSliceNormalVector(Header)
%==========================================================================
function nrm = ReadSliceNormalVector(Header)
str = Header.CSAImageHeaderInfo;
val = GetNumaris4Val(str,'SliceNormalVector');
for i=1:3
    nrm(i,1) = sscanf(val(i,:),'%g');
end


%==========================================================================
% function n = ReadNumberOfImagesInMosaic(Header)
%==========================================================================
function n = ReadNumberOfImagesInMosaic(Header)
str = Header.CSAImageHeaderInfo;
val = GetNumaris4Val(str,'NumberOfImagesInMosaic');
n   = sscanf(val','%d');
if isempty(n), n = []; end


%==========================================================================
% function dim = ReadAcquisitionMatrixText(Header)
%==========================================================================
function dim = ReadAcquisitionMatrixText(Header)
str = Header.CSAImageHeaderInfo;
val = GetNumaris4Val(str,'AcquisitionMatrixText');
dim = sscanf(val','%d*%d')';
if length(dim)==1
    dim = sscanf(val','%dp*%d')';
end
if isempty(dim), dim=[]; end


%==========================================================================
% function val = GetNumaris4Val(str,name)
%==========================================================================
function val = GetNumaris4Val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str)
    if strcmp(deblank(str(i).name),name)
        for j=1:str(i).nitems
            if  str(i).item(j).xx(1)
                val = [val {str(i).item(j).val}];
            end
        end
        break;
    end
end
val = strvcat(val{:});


%==========================================================================
% function val = GetNumaris4NumVal(str,name)
%==========================================================================
function val = GetNumaris4NumVal(str,name)
val1 = GetNumaris4Val(str,name);
val  = zeros(size(val1,1),1);
for k = 1:size(val1,1)
    val(k)=str2num(val1(k,:));
end


%==========================================================================
% function fname = getfilelocation(Header,RootDirectory,prefix,format,OutputDirectory)
%==========================================================================
function fname = getfilelocation(Header,RootDirectory,prefix,format,OutputDirectory)

if nargin < 3
    prefix = 'f';
end

if strncmp(RootDirectory,'ice',3)
    RootDirectory = RootDirectory(4:end);
    imtype = textscan(Header.ImageType,'%s','delimiter','\\');
    try
        imtype = imtype{1}{3};
    catch
        imtype = '';
    end
    prefix = [prefix imtype GetNumaris4Val(Header.CSAImageHeaderInfo,'ICE_Dims')];
end

if isfield(Header,'PatientID'),         PatientID         = deblank(Header.PatientID); else PatientID         = 'anon'; end
if isfield(Header,'EchoNumbers'),       EchoNumbers       = Header.EchoNumbers;        else EchoNumbers       = 0;      end
if isfield(Header,'SeriesNumber'),      SeriesNumber      = Header.SeriesNumber;       else SeriesNumber      = 0;      end
if isfield(Header,'AcquisitionNumber'), AcquisitionNumber = Header.AcquisitionNumber;  else AcquisitionNumber = 0;      end
if isfield(Header,'InstanceNumber'),    InstanceNumber    = Header.InstanceNumber;     else InstanceNumber    = 0;      end

ImTyp = '';
if isfield(Header,'GE_ImageType')
    if numel(Header.GE_ImageType)==1
        switch Header.GE_ImageType
        case 1
            ImTyp = '-Phase';
        case 2
            ImTyp = '-Real';
        case 3
            ImTyp = '-Imag';
        end
    end
end

% To use ICE Dims systematically in file names in order to avoid
% overwriting uncombined coil images, which have identical file name
% otherwise)
try 
    ICE_Dims = GetNumaris4Val(Header.CSAImageHeaderInfo,'ICE_Dims');
    % extract ICE dims as an array of numbers (replace 'X' which is for
    % combined images by '-1' first): 
    CHA = sscanf(strrep(ICE_Dims,'X','-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
    if CHA(1)>0
        CHA = sprintf('%.3d',CHA(1));
    else 
        CHA = '';
    end
catch
    CHA = '';
end

if strcmp(RootDirectory, 'flat')
    % Standard SPM file conversion
    %----------------------------------------------------------------------
    if CheckFields(Header,'SeriesNumber','AcquisitionNumber')
        if CheckFields(Header,'EchoNumbers')
            if ~isempty(CHA)
                fname = sprintf('%s%s-%.4d-%.5d-%.6d-%.2d-%s%s.%s', prefix, StripUnwantedChars(PatientID),...
                                SeriesNumber, AcquisitionNumber, InstanceNumber, EchoNumbers, CHA, ImTyp, format);
            else
                fname = sprintf('%s%s-%.4d-%.5d-%.6d-%.2d%s.%s', prefix, StripUnwantedChars(PatientID),...
                                SeriesNumber, AcquisitionNumber, InstanceNumber, EchoNumbers, ImTyp, format);
            end
        else
            if ~isempty(CHA)
                fname = sprintf('%s%s-%.4d-%.5d-%.6d-%s%s.%s', prefix, StripUnwantedChars(PatientID),...
                                SeriesNumber, AcquisitionNumber, InstanceNumber, CHA, ImTyp, format);
            else
                fname = sprintf('%s%s-%.4d-%.5d-%.6d%s.%s', prefix, StripUnwantedChars(PatientID),...
                                SeriesNumber, AcquisitionNumber, InstanceNumber, ImTyp, format);
            end
        end
    else
        fname = sprintf('%s%s-%.6d%s.%s',prefix, ...
                        StripUnwantedChars(PatientID), InstanceNumber, ImTyp, format);
    end

    fname = fullfile(OutputDirectory,fname);
    return;
end

% more fancy stuff - sort images into subdirectories
if isfield(Header,'StudyTime')
    m = sprintf('%02d', floor(rem(Header.StudyTime/60,60)));
    h = sprintf('%02d', floor(Header.StudyTime/3600));
else
    m = '00';
    h = '00';
end
if isfield(Header,'AcquisitionTime'),   AcquisitionTime   = Header.AcquisitionTime;            else AcquisitionTime   = 100;       end
if isfield(Header,'StudyDate'),         StudyDate         = Header.StudyDate;                  else StudyDate         = 100;       end % Obscure Easter Egg
if isfield(Header,'SeriesDescription'), SeriesDescription = deblank(Header.SeriesDescription); else SeriesDescription = 'unknown'; end
if isfield(Header,'ProtocolName')
    ProtocolName = deblank(Header.ProtocolName);
else
    if isfield(Header,'SequenceName')
        ProtocolName = deblank(Header.SequenceName);
    else
        ProtocolName='unknown';
    end
end

studydate = sprintf('%s_%s-%s', datestr(StudyDate,'yyyy-mm-dd'), h,m);
switch RootDirectory
    case {'date_time','series'}
        id = studydate;
    case {'patid', 'patid_date', 'patname'}
        id = StripUnwantedChars(PatientID);
end
serdes   = strrep(StripUnwantedChars(SeriesDescription), StripUnwantedChars(ProtocolName),'');
protname = sprintf('%s%s_%.4d', StripUnwantedChars(ProtocolName), serdes, SeriesNumber);
switch RootDirectory
    case 'date_time'
        dname = fullfile(OutputDirectory, id, protname);
    case 'patid'
        dname = fullfile(OutputDirectory, id, protname);
    case 'patid_date'
        dname = fullfile(OutputDirectory, id, studydate, protname);
    case 'series'
        dname = fullfile(OutputDirectory, protname);
    otherwise
        error('unknown file root specification');
end
if ~exist(dname,'dir')
    mkdir_rec(dname);
end

% some non-product sequences on SIEMENS scanners seem to have problems
% with image numbering in MOSAICs - doublettes, unreliable ordering
% etc. To distinguish, always include Acquisition time in image name
sa    = sprintf('%02d', floor(rem(AcquisitionTime,60)));
ma    = sprintf('%02d', floor(rem(AcquisitionTime/60,60)));
ha    = sprintf('%02d', floor(AcquisitionTime/3600));
fname = sprintf('%s%s-%s%s%s-%.5d-%.5d-%d%s.%s', prefix, id, ha, ma, sa, ...
        AcquisitionNumber, InstanceNumber, EchoNumbers, ImTyp, format);
fname = fullfile(dname, fname);


%==========================================================================
% function suc = mkdir_rec(str)
%==========================================================================
function suc = mkdir_rec(str)
% works on full pathnames only
if str(end) ~= filesep, str = [str filesep]; end
pos = strfind(str,filesep);
suc = zeros(1,length(pos));
for g=2:length(pos)
    if ~exist(str(1:pos(g)-1),'dir')
        suc(g) = mkdir(str(1:pos(g-1)-1),str(pos(g-1)+1:pos(g)-1));
    end
end


%==========================================================================
% function ret = ReadAscconv(Header)
%==========================================================================
function ret = ReadAscconv(Header)
% In SIEMENS data, there is an ASCII text section with
% additional information items. This section starts with a code
% ### ASCCONV BEGIN <some version string> ###
% and ends with
% ### ASCCONV END ###
% It is read by spm_dicom_headers into an entry 'MrProtocol' in
% CSASeriesHeaderInfo or into an entry 'MrPhoenixProtocol' in
% CSAMiscProtocolHeaderInfoVA or CSAMiscProtocolHeaderVB.
% The additional items are assignments in C syntax, here they are just
% translated according to
% [] -> ()
% "  -> '
% 0xX -> hex2dec('X')
% and collected in a struct.
% In addition, there seems to be "__attribute__" meta information for some
% items. All "field names" starting with "_" are currently silently ignored.
ret=struct;

% get ascconv data
if isfield(Header, 'CSAMiscProtocolHeaderInfoVA')
    X = GetNumaris4Val(Header.CSAMiscProtocolHeaderInfoVA,'MrProtocol');
elseif isfield(Header, 'CSAMiscProtocolHeaderInfoVB')
    X = GetNumaris4Val(Header.CSAMiscProtocolHeaderInfoVB,'MrPhoenixProtocol');
elseif isfield(Header, 'CSASeriesHeaderInfo')
    X = GetNumaris4Val(Header.CSASeriesHeaderInfo,'MrProtocol');
else
    return;
end

X = regexprep(X,'^.*### ASCCONV BEGIN [^#]*###(.*)### ASCCONV END ###.*$','$1');

if ~isempty(X)
    tokens = textscan(char(X),'%s', ...
        'delimiter',char(10));
    tokens{1}=regexprep(tokens{1},{'\[([0-9]*)\]','"(.*)"','^([^"]*)0x([0-9a-fA-F]*)','#.*','^.*\._.*$'},{'($1+1)','''$1''','$1hex2dec(''$2'')','',''});
    % If everything would evaluate correctly, we could use
    % eval(sprintf('ret.%s;\n',tokens{1}{:}));
    for k = 1:numel(tokens{1})
        if ~isempty(tokens{1}{k})
            try
                [tlhrh] = regexp(tokens{1}{k}, '(?:=)+', 'split', 'match');
                [tlh]   = regexp(tlhrh{1}, '(?:\.)+', 'split', 'match');
                tlh = cellfun(@genvarname, tlh, 'UniformOutput',false);
                tlh = sprintf('.%s', tlh{:});
                eval(sprintf('ret%s = %s;', tlh, tlhrh{2}));
            catch
                disp(['AscConv: Error evaluating ''ret.' tokens{1}{k} ''';']);
            end
        end
    end
end


%==========================================================================
% function Datatype = DetermineDatatype(Header)
%==========================================================================
function Datatype = DetermineDatatype(Header)
% Determine what datatype to use for NIfTI images
BigEndian = spm_platform('bigend');
if Header.BitsStored>16
    if Header.PixelRepresentation
        Datatype = [spm_type( 'int32') BigEndian];
    else
        Datatype = [spm_type('uint32') BigEndian];
    end
else
    if Header.PixelRepresentation 
        Datatype = [spm_type( 'int16') BigEndian];
    else
        Datatype = [spm_type('uint16') BigEndian];
    end
end


%==========================================================================
% function fspe = ConvertMultiframes(Headers, RootDirectory, format, OutputDirectory, meta)
%==========================================================================
function fspe = ConvertMultiframes(Headers, RootDirectory, format, OutputDirectory, meta)
fspe = {};
DicomDictionary = ReadDicomDictionary;
for i=1:numel(Headers)
    out  = ConvertMultiframe(Headers{i}, DicomDictionary, RootDirectory, format,OutputDirectory,meta);
    fspe = [fspe(:); out(:)];
end


%==========================================================================
% function out = ConvertMultiframe(Header, DicomDictionary, RootDirectory, format, OutputDirectory, meta)
%==========================================================================
function out = ConvertMultiframe(Header, DicomDictionary, RootDirectory, format, OutputDirectory, meta)
out      = {};
DimOrg   = GetDimOrg(Header, DicomDictionary);
FGS      = GetFGS(Header, DimOrg);
N        = numel(FGS);
fname0   = getfilelocation(Header, RootDirectory,'MF',format,OutputDirectory);
[pth,nam,ext0] = fileparts(fname0);
fname0   = fullfile(pth,nam);

if ~isfield(FGS,'ImageOrientationPatient') || ~isfield(FGS,'PixelSpacing') || ~isfield(FGS,'ImagePositionPatient')
    fprintf('"%s" does not seem to have positional information.\n', Header.Filename);
    return;
end

%if isfield(FGS, 'DimensionIndexValues')
%    disp(cat(1, FGS.DimensionIndexValues))
%end

% Read the actual image data
volume = ReadImageData(Header);

% Sort into unique files according to image orientations etc
stuff = [cat(2, FGS.ImageOrientationPatient); cat(2, FGS.PixelSpacing)]';
if isfield(FGS, 'StackID')
    stuff = [stuff double(cat(1, FGS.StackID))];
end
ds    = find(any(diff(stuff,1,1)~=0,2));
ord   = sparse([],[],[],N,numel(ds)+1);
start = 1;
for i=1:numel(ds)
    ord(start:ds(i),i) = 1;
    start = ds(i)+1;
end
ord(start:size(stuff,1),numel(ds)+1) = 1;
for n=1:size(ord,2)
    % Identify all slices that should go into the same output file
    %ind  = find(all(bsxfun(@eq,stuff,u(n,:)),2));
    ind = find(ord(:,n));
    this = FGS(ind); 
    if size(ord,2)>1
        fname = sprintf('%s_%-.3d',fname0,n);
    else
        fname = fname0;
    end


    % Orientation information
    %----------------------------------------------------------------------
    % Axial Analyze voxel co-ordinate system:
    % x increases     right to left
    % y increases posterior to anterior
    % z increases  inferior to superior

    % DICOM patient co-ordinate system:
    % x increases     right to left
    % y increases  anterior to posterior
    % z increases  inferior to superior

    % T&T co-ordinate system:
    % x increases      left to right
    % y increases posterior to anterior
    % z increases  inferior to superior


    if any(sum(diff(cat(2,this.ImageOrientationPatient),1,2).^2,1)>0.001)
        fprintf('"%s" contains irregularly oriented slices.\n', Header.Filename);
        % break % Option to skip writing out the NIfTI image
    end
    ImageOrientationPatient = this(1).ImageOrientationPatient(:);

    if any(sum(diff(cat(2,this.PixelSpacing),1,2).^2,1)>0.001)
        fprintf('"%s" contains slices with irregularly spaced pixels.\n', Header.Filename);
        % break % Option to skip writing out the NIfTI image
    end
    PixelSpacing            = this(1).PixelSpacing(:);

    R   = [reshape(ImageOrientationPatient,3,2)*diag(PixelSpacing); 0 0];

    % Determine the order of the slices by sorting their positions according to where they project
    % on a vector orthogonal to the image plane
    orient      = reshape(this(1).ImageOrientationPatient,[3 2]);
    orient(:,3) = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
        
    Positions  = cat(2,this.ImagePositionPatient)';
    [proj,slice_order,inv_slice_order] = unique(round(Positions*orient(:,3)*100)/100);
    if any(abs(diff(diff(proj),1,1))>0.025)
        problem1 = true;
    else
        problem1 = false;
    end
 
    inv_time_order = ones(numel(this),1);
    for i=1:numel(slice_order)
        ind = find(inv_slice_order == i);
        if numel(ind)>1
            if isfield(this,'TemporalPositionIndex')
                sort_on = cat(2,this(ind).TemporalPositionIndex)';
            else
                sort_on = ind;
            end
            [unused, tsort]     = sort(sort_on);
            inv_time_order(ind) = tsort';
        end
    end

    
    % Image dimensions
    %----------------------------------------------------------------------
    nc   = Header.Columns;
    nr   = Header.Rows;
    dim  = [nc nr 1 1];
    dim(3) = max(inv_slice_order);
    dim(4) = max(inv_time_order);

    problem2 = false;
    for i=1:max(inv_time_order)
        if sum(inv_time_order==i)~=dim(3)
            problem2 = true;
        end
    end
    if problem1 || problem2
        fname = [fname '-problem'];
    end
    if problem1
        fprintf('"%s" contains irregularly spaced slices with the following spacings:\n', Header.Filename);
        spaces = unique(round(diff(proj)*100)/100);
        fprintf(' %g', spaces);
        fprintf('\nSee %s%s\n\n',fname,ext0);
        % break % Option to skip writing out the NIfTI image
    end
    if problem2
        fprintf('"%s" has slices missing.\nSee the result in %s%s\n\n', Header.Filename,fname,ext0);
        % break % Option to skip writing out the NIfTI image
    end
 
    if dim(3)>1
        y1 = [this(slice_order(  1)).ImagePositionPatient(:); 1];
        y2 = [this(slice_order(end)).ImagePositionPatient(:); 1];
        x2 = [1; 1; dim(3); 1];
    else
        orient      = reshape(ImageOrientationPatient,[3 2]);
        orient(:,3) = null(orient');
        if det(orient)<0, orient(:,3) = -orient(:,3); end
        if isfield(Header,'SliceThickness'), z = Header.SliceThickness; else z = 1; end
        y1 = [this(1).ImagePositionPatient(:); 1];
        y2 = [orient*[0;0;z];0];
        x2 = [0;0;1;0];
    end

    x1             = [1;1;1;1];
    DicomToPatient = [y1 y2 R]/[x1 x2 eye(4,2)];
    AnalyzeToDicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
    PatientToTal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions
    mat            = PatientToTal*DicomToPatient*AnalyzeToDicom;
    flip_lr        = det(mat(1:3,1:3))>0;

    if flip_lr
        mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
    end

    % Possibly useful information for descrip field
    %----------------------------------------------------------------------
    if isfield(Header,'AcquisitionTime')
        tim = datevec(Header.AcquisitionTime/(24*60*60));
    elseif isfield(Header,'StudyTime')
        tim = datevec(Header.StudyTime/(24*60*60));
    elseif isfield(Header,'ContentTime') 
        tim = datevec(Header.ContentTime/(24*60*60));      
    else
        tim = '';
    end
    if ~isempty(tim), tim = sprintf(' %d:%d:%.5g', tim(4),tim(5),tim(6)); end

    if isfield(Header,'AcquisitionDate') 
        day = datestr(Header.AcquisitionDate);
    elseif isfield(Header,'StudyDate')
        day = datestr(Header.StudyDate);
    elseif isfield(Header,'ContentDate')
        day = datestr(Header.ContentDate);
    else
        day = '';
    end
    when = [day tim]; 

    if CheckFields(Header,'MagneticFieldStrength','MRAcquisitionType',...
                   'ScanningSequence','RepetitionTime','EchoTime','FlipAngle')
        if isfield(Header,'ScanOptions')
            ScanOptions = Header.ScanOptions;
        else
            ScanOptions = 'no';
        end
        modality = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg/SO=%s',...
                           Header.MagneticFieldStrength, Header.MRAcquisitionType,...
                           deblank(Header.ScanningSequence),...
                           Header.RepetitionTime, Header.EchoTime, Header.FlipAngle,...
                           ScanOptions);
    else
         modality = [Header.Modality ' ' Header.ImageType];
    end
    descrip = [modality ' ' when];



    % Sort out datatype, as well as any scalefactors or intercepts
    %----------------------------------------------------------------------
    pinfos = [ones(length(this),1) zeros(length(this),1)];
    for i=1:length(this)
        if isfield(this(i),'RescaleSlope'),     pinfos(i,1) = this(i).RescaleSlope;     end
        if isfield(this(i),'RescaleIntercept'), pinfos(i,2) = this(i).RescaleIntercept; end

        % Philips do things differently. The following is for using their scales instead.
        %     Chenevert, Thomas L., et al. "Errors in quantitative image analysis due to
        %     platform-dependent image scaling." Translational oncology 7.1 (2014): 65-71.
        if isfield(this(i),'MRScaleSlope'),     pinfos(i,1) = 1/this(i).MRScaleSlope;                 end
        if isfield(this(i),'MRScaleIntercept'), pinfos(i,2) =  -this(i).MRScaleIntercept*pinfos(i,1); end

    end
    
    % The following requires more testing. Only tested on the following 
    % multiframe datasets (all Philips MRI data):
    % - EPI series (either GRE-EPI (fMRI) or SE-EPI (DWI)),
    % - structural mono-volume images (FLAIR, T1FFE),
    % - structural multi-echo images (MPM protocol)
    % - B0 mapping (two volumes including magnitude and phase in
    %   milliradians).
    % The first three cases above correspond to 
    % any(any(diff(pinfos,1))) == false (identical scaling for all frames,
    % no problem there), while the latter corresponds to
    % any(any(diff(pinfos,1))) == true (variable scaling and problems).  
    if ~any(any(diff(pinfos,1)))
        % Same slopes and intercepts for all slices
        dt    = DetermineDatatype(Header);
        pinfo = pinfos(1,:);
    else
        % Variable slopes and intercept (maybe PET/SPECT) - or
        % magnitude/phase data, or...?

        % If variable slopes and intercepts in the multiframe dataset, 
        % each volume is saved as a separate NIfTI file with appropriate
        % slope and intercept. NB: we assume that the same pinfo is used
        % for all frames in a given volume - might not be the case?
        dt = DetermineDatatype(Header); 
        dt    = 'int16-be';
        pinfo = [];
        for cv = 1:max(inv_time_order)
            ind = find(inv_time_order==cv);
            pinfo = [pinfo; pinfos(ind(1),:)]; 
        end
        
        % Ensure random numbers are reproducible (see later)
        % when intensities are dithered to prevent aliasing effects.
        rand('state',0);
    end

    % Define output NIfTI files
    %----------------------------------------------------------------------
    fname0 = fname;
    for cNii = 1:size(pinfo,1)
        % Define output file name
        if size(pinfo,1)>1
            fname = sprintf('%s-%.4d',fname0,cNii);
            dim = dim(1:3);
        end
        
        % Create nifti file
        Nii{cNii}      = nifti; %#ok<*AGROW>
        fname          = sprintf('%s%s', fname, ext0);
        Nii{cNii}.dat  = file_array(fname, dim, dt, 0, pinfo(cNii,1), pinfo(cNii,2));
        Nii{cNii}.mat  = mat;
        Nii{cNii}.mat0 = mat;
        Nii{cNii}.mat_intent  = 'Scanner';
        Nii{cNii}.mat0_intent = 'Scanner';
        Nii{cNii}.descrip     = descrip;
        create(Nii{cNii});
    
        % for json metadata, need to sort out the portion of the Header
        % that corresponds to the current Nii volume (if multiscale &
        % multiframe data -> 3D volumes): 
        if meta
            cHeader{cNii} = Header;
            if size(pinfo,1)>1
                cHeader{cNii} = rmfield(cHeader{cNii},'PerFrameFunctionalGroupsSequence');
                ind = find(inv_time_order==cNii);
                cHeader{cNii}.NumberOfFrames = length(ind);
                for cF = 1:length(ind)
                    cHeader{cNii}.PerFrameFunctionalGroupsSequence{cF} = Header.PerFrameFunctionalGroupsSequence{ind(cF)};
                end
            end
            Nii{cNii} = spm_dicom_metadata(Nii{cNii},cHeader{cNii});
        end
        
        Nii{cNii}.dat(end,end,end,end,end) = 0;
    end

    % Write the image volume
    %----------------------------------------------------------------------
    spm_progress_bar('Init',length(this),['Writing ' fname], 'Planes written');

    for i=1:length(this)

        plane = volume(:,:,this(i).number);

        if any(any(diff(pinfos,1)))
            % This is to prevent aliasing effects in any subsequent histograms
            % of the data (eg for mutual information coregistration).
            % It's a bit inelegant, but probably necessary for when slices are
            % individually rescaled.
            plane = double(plane) + rand(size(plane)) - 0.5;
        end

        if pinfos(i,1)~=1, plane = plane*pinfos(i,1); end
        if pinfos(i,2)~=0, plane = plane+pinfos(i,2); end

        plane = fliplr(plane);
        if flip_lr, plane = flipud(plane); end

        z = inv_slice_order(i);
        t = inv_time_order(i);
        if size(pinfo,1)>1
            Nii{t}.dat(:,:,z) = plane;
        else
            Nii{1}.dat(:,:,z,t) = plane;
        end
        spm_progress_bar('Set',i);
    end

    for cNii = 1:size(pinfo,1)
        out = [out; {Nii{cNii}.dat.fname}];
    end
    spm_progress_bar('Clear');
end


%==========================================================================
% function DimOrg = GetDimOrg(Header, DicomDictionary)
%==========================================================================
function DimOrg = GetDimOrg(Header, DicomDictionary)
DimOrg = struct('DimensionIndexPointer',{},'FunctionalGroupPointer',{});

if isfield(Header,'DimensionIndexSequence')
    for i=1:numel(Header.DimensionIndexSequence)
        DIP = Header.DimensionIndexSequence{i}.DimensionIndexPointer;
        ind = find(DicomDictionary.group==DIP(1) & DicomDictionary.element == DIP(2));
        if numel(ind)==1
            DimOrg(i).DimensionIndexPointer = DicomDictionary.values(ind).name;
        else
            if rem(DIP(1),2)
                DimOrg(i).DimensionIndexPointer = sprintf('Private_%.4x_%.4x',DIP);
            else
                DimOrg(i).DimensionIndexPointer = sprintf('Tag_%.4x_%.4x',DIP);
            end
        end

        FGP = Header.DimensionIndexSequence{i}.FunctionalGroupPointer;
        ind = find(DicomDictionary.group==FGP(1) & DicomDictionary.element == FGP(2));
        if numel(ind)==1
            DimOrg(i).FunctionalGroupPointer     = DicomDictionary.values(ind).name;
        else
            if rem(DIP(1),2)
                DimOrg(i).FunctionalGroupPointer = sprintf('Private_%.4x_%.4x',DIP);
            else
                DimOrg(i).FunctionalGroupPointer = sprintf('Tag_%.4x_%.4x',DIP);
            end
        end
    end
end
return


%==========================================================================
% FGS = GetFGS(Header, DimOrg)
%==========================================================================
function FGS = GetFGS(Header, DimOrg)
FGS = struct;
if isfield(Header,'PerFrameFunctionalGroupsSequence')
    % For documentation, see C.7.6.16 "Multi-frame Functional Groups Module"
    % of DICOM Standard PS 3.3 - 2003, Page 261.

    % Multiframe DICOM
    N = numel(Header.PerFrameFunctionalGroupsSequence);
    if isfield(Header,'NumberOfFrames') && Header.NumberOfFrames ~= N
        fprintf('"%s" has incompatible numbers of frames.', FGS.Filename);
    end

    % List of fields and subfields of those fields to retrieve information from.
    Macros = {'PixelMeasuresSequence',{'PixelSpacing','SliceThickness'}
              'FrameContentSequence',{'Frame Acquisition Number','FrameReferenceDatetime',...
                                      'FrameAcquisitionDatetime','FrameAcquisitionDuration',...
                                      'CardiacCyclePosition','RespiratoryCyclePosition',...
                                      'DimensionIndexValues','TemporalPositionIndex','Stack ID',...
                                      'InStackPositionNumber','FrameComments'}
              'PlanePositionSequence',{'ImagePositionPatient'}
              'PlaneOrientationSequence',{'ImageOrientationPatient'}
              'ReferencedImageSequence',{'ReferencedSOPClassUID','ReferencedSOPInstanceUID',...
                                         'ReferencedFrameNumber','PurposeOfReferenceCode'}
              'PixelValueTransformationSequence',{'RescaleIntercept','RescaleSlope','RescaleType'}
              % for specific PHILIPS rescaling
              'PhilipsSequence_2005_140f',{'MRScaleSlope','MRScaleIntercept'}};

    if isfield(Header,'SharedFunctionalGroupsSequence')
        for d=1:numel(DimOrg)
            if isfield(Header.SharedFunctionalGroupsSequence{1}, DimOrg(d).FunctionalGroupPointer) && ...
               isfield(Header.SharedFunctionalGroupsSequence{1}.(DimOrg(d).FunctionalGroupPointer){1}, DimOrg(d).DimensionIndexPointer)
                for n=N:-1:1
                    FGS(n).(DimOrg(d).DimensionIndexPointer) = Header.SharedFunctionalGroupsSequence{1}.(DimOrg(d).FunctionalGroupPointer){1}.(DimOrg(d).DimensionIndexPointer);
                end
            end
        end

        for k1=1:size(Macros,1)
            if isfield(Header.SharedFunctionalGroupsSequence{1}, Macros{k1,1})
                for k2=1:numel(Macros{k1,2})
                    if (numel(Header.SharedFunctionalGroupsSequence{1}.(Macros{k1,1}))>0) &&...
                            isfield(Header.SharedFunctionalGroupsSequence{1}.(Macros{k1,1}){1}, Macros{k1,2}{k2})
                        for n=N:-1:1
                            FGS(n).(Macros{k1,2}{k2}) = Header.SharedFunctionalGroupsSequence{1}.(Macros{k1,1}){1}.(Macros{k1,2}{k2});
                        end
                    end
                end
            end
        end
    end

    for n=1:N
        FGS(n).number = n;
        F = Header.PerFrameFunctionalGroupsSequence{n};
        for d=1:numel(DimOrg)
            if isfield(F,DimOrg(d).FunctionalGroupPointer) && ...
               isfield(F.(DimOrg(d).FunctionalGroupPointer){1}, DimOrg(d).DimensionIndexPointer)
                FGS(n).(DimOrg(d).DimensionIndexPointer) = F.(DimOrg(d).FunctionalGroupPointer){1}.(DimOrg(d).DimensionIndexPointer);
            end
        end

        for k1=1:size(Macros,1)
            if isfield(F,Macros{k1,1}) && (numel(F.(Macros{k1,1}))>0)
                for k2=1:numel(Macros{k1,2})
                    if isfield(F.(Macros{k1,1}){1},Macros{k1,2}{k2})
                        FGS(n).(Macros{k1,2}{k2}) = F.(Macros{k1,1}){1}.(Macros{k1,2}{k2});
                    end
                end
            end
        end
    end
else
    error('"%s" is not multiframe.', Header.Filename);
end


%==========================================================================
% function DicomDictionary = ReadDicomDictionary(DictionaryFile)
%==========================================================================
function DicomDictionary = ReadDicomDictionary(DictionaryFile)
if nargin<1, DictionaryFile = 'spm_dicom_dict.mat'; end
try
    DicomDictionary = load(DictionaryFile);
catch Problem
    fprintf('\nUnable to load the file "%s".\n', DictionaryFile);
    rethrow(Problem);
end

