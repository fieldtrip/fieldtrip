function [MRItag,MRIdata]=readCTFMRI(fileName);

%  Version 1.2: 25 April 2007   Module readCPersist changed and removed from this listing.
%  Version 1.1: 19 April 2007   No changes since v1.0

%  Reads a CTFMRI v4 file.
%  For file format and a defineition fo CPersist objects, see document 
%  "CTF MEG FIle Formats", PN900-0088

%  Input: fileName : Name of MRI file  including path and extension.

%  Outputs: MRItags : CPersist tags describing the MRI.
%                     The start and stop tags have been removed.
%                        name='WS1_', type=0, data=[]
%                        name='EndOfParmaters', type -1, data=[]

%           MRIdata: 256x256x256: int16 array.
%               MRIdata(j,k,n) = pixel(j,k) of slice n.
%                               j = coronal slice.  j=1:anterior, j=256:posterior
%                               k = axial slice.  k=1:superior, k=256:inferior
%                               n = sagittal slice.  n=1:left, n=256:right

%  Calls function readCPersist (not included in this listing).

persistent printWarning


if nargin==0 & nargout==0
  fprintf(['readCTFMRI:  Version 1.2   25 April 2007   Reads MRIs in CTFMRI v4 format.\n',...
      '\t[MRItag,MRIdata]=readCTFMRI(MRIfilename);\n',...
      '\t\tMRIfilename = complete name including path and extension .mri\n',...
      '\t\tMRItag= list of File descriptors.  See document "CTF MEG FIle Formats", PN900-0088\n',...
      '\t\tMRIdata = MRI data array.  Precision=int16 (2 bytes per pixel)\n\n']);
  return
end

MRItag=struct([]);  % Force failure in calling routines
MRIdata=int16([]);

if nargin~=1 | nargout~=2
  fprintf(['readCTFMRI: Found %d input arguments and %d outputs. ',...
      'Must have 1 inputs (fileName) and 2 outputs [MRItags,MRIdata].\n'],nargin,nargout);
  return
elseif ~ischar(fileName) | isempty(fileName) 
  fprintf('readCTFMRI : Filename is not character or is empty.\n');
  whos fileName;
  return
elseif exist(fileName)~=2
  fprintf('readCTFMRI: Cannot find file %s.\n',fileName);
  return
end

MRItag=readCPersist(fileName);
if isempty(MRItag)
  fprintf('readCTFMRI:  File %s\n\t\t\t cannot be read by readCPersist.\n',fileName);
  fprintf('\t\t\tIs it in CPersist format?\n\n');
  MRItag=struct([]);
  MRIdata=int16([]);
  return
end

%  DO some checks to verify that this really is a CTF MRI.
name=char([]);
for k=1:length(MRItag);
  name=strvcat(name,MRItag(k).name);
end
isCTFMRI=(~isempty(strmatch('_CTFMRI_',name)));
if isCTFMRI
  sizeTag=strmatch('_CTFMRI_SIZE',name,'exact');
  if length(sizeTag)~=1
    isCTFMRI=0;
  else
    MRIsize=MRItag(sizeTag).data;
    nSlice=MRIsize;
  end
  clear sizeTag;
end
if isCTFMRI
  sliceTag=strmatch('_CTFMRI_SLICE_DATA#',name)';
  isCTFMRI=(length(sliceTag)==nSlice);
end
if ~isCTFMRI
  fprintf(['readCTFMRI:  File %s is a CPersist file\n',...
      '\t\t\tbut it is not in CTFMRI format.\n'],fileName);
  MRIdata=int16([]);
  return
end
clear isCTFMRI;

if isempty(printWarning)
  fprintf(['\nreadCTFMRI: You are reading a CTF-format MRI for use with a software-application\n',...
      '\ttool that is not manufactured by VSM MedTech Ltd. and has not received marketing\n',...
      '\tclearance for clinical applications.  If the MRI is processed by this tool,\n',...
      '\tit should not be later employed for clinical and/or diagnostic purposes.\n']);
  printWarning=1;
end

MRIdata=int16(zeros(MRIsize,MRIsize,1));
for k=1:nSlice
  MRIdata(:,:,k)=reshape(int16(MRItag(sliceTag(k)).data),MRIsize,MRIsize);
end
% Keep only the useful tag information.
MRItag=MRItag(2:(sliceTag(1)-1));
return
%%%%%%%%%%%%%%  end of readCTFMRI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%