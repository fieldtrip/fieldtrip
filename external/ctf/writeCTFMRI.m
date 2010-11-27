function writeCTFMRI(fileName,MRItags,MRIdata);

% Version 1.2   25 April 2007   Module writeCPersist changed, and removed from this text
%                               file.
% Version 1.1   19 April 2007 : No changes from v1.0
% Version 1.0   27 Oct. 2006

%  Write a CTF MRI in v4 format.
%  For file format and a defineition fo CPersist objects, see document 
%  "CTF MEG FIle Formats", PN900-0088

%  Input: fileName : Name of MRI file  including path and extension.

%  Outputs: MRItags : CPersist tags describing the MRI in the format of the structure array
%                     produced by readMRI.  There must be no start and stop tags.
%                     Clnical use messages are added to _PATIENT_ID, _PATIENT_NAME,
%                     _STUDY_COMMENTS, _CTFMRI_COMMENTS.
%                      I.e  name='WS1_', type=0, data=[] and
%                           name='EndOfParameters', type -1, data=[]

%           MRIdata: 256x256x256 array.  No value may exceed 2^15-1 since the data are
%                                        converted to int16.
%               MRIdata(j,k,n) = pixel(j,k) of slice n.
%                               j = coronal slice.  j=1:anterior, j=256:posterior
%                               k = axial slice.  k=1:superior, k=256:inferior
%                               n = sagittal slice. j=1:left, j=256:right

%  Calls: writeCPersist (not included in this listing)
%         update_MRI_descriptors (included in this listing)

persistent printWarning

clinical_use_message='NOT FOR CLINICAL USE';
creatorSoftware='writeCTFMRI.m';  % Added to .infods file

NpxNot256Warning=1;  % Print a warning if the no. of slices is not 256

if nargin==0 & nargout==0
  fprintf(['writeCTFMRI:  Version 1.2   25 April 2007   Creates MRIs in CTFMRI v4 format.\n',...
      '\twriteCTFMRI(MRIfilename,MRItags,MRIdata);\n',...
      '\t\tMRIfilename = complete name including path and extension .mri\n',...
      '\t\tMRItags= structure array listing descriptors in the format produced',...
      'by function readMRI.\n',...
      '\t\tMRIdata = MRI data array.\n\n',...
      '\tSee document "CTF MEG FIle Formats", PN900-0088\n']);
  return
end

if nargin~=3
  fprintf(['writeCTFMRI: Found %d input arguments.  ',...
      'Must have 3 inputs (filename,MRItags,MRIdata).\n'],nargin);
  return
elseif ~ischar(fileName) | isempty(fileName) | ~isstruct(MRItags) | isempty(MRItags) |...
    ~isnumeric(MRIdata) | isempty(MRIdata)
  fprintf('writeCTFMRI : Wrong argument type, or argument(s) empty.\n');
  whos fileName MRItags MRIdata;
  return
elseif ~isfield(MRItags,'name') | ~isfield(MRItags,'type') | ~isfield(MRItags,'data')
  fprintf('writeCTFMRI: MRItags does not have the correct fields (name, type,data).\n');
  MRItags
  return
elseif ndims(MRIdata)~=3 | ~all(size(MRIdata)==size(MRIdata,1))
  fprintf('writeCTFMRI: size(MRIdata)=[');fprintf(' %d',size(MRIdata));...
      fprintf(']   Must be Npx*[1 1 1] in CTF MRI format.\n');
  return
elseif exist(fileName)==2
  fprintf('writeCTFMRI: File %s already exists.\n',fileName);
  return
end

%  Verify that this really is a CTF MRI.
name=char([]);
for k=1:length(MRItags);
  name=strvcat(name,MRItags(k).name);
end
isCTFMRI=(~isempty(strmatch('_CTFMRI_',name)));

MRIsize=size(MRIdata,1);  % assume square slices
nSlice=size(MRIdata,3);
if NpxNot256Warning & MRIsize~=256
  fprintf(['\nwriteCTFMRI: size(MRIdata)=%d*[1 1 1].  CTF MRI format standard is ',...
      '256*[1 1 1].\n\t\t\tThis file may not work with MRIViewer.\n\n'],MRIsize);
end
if isCTFMRI
  sizeTag=strmatch('_CTFMRI_SIZE',name,'exact');
  if length(sizeTag)~=1
    isCTFMRI=0;
  else
    MRItags(sizeTag).data=MRIsize;
  end
  clear sizeTag;
end
if ~isCTFMRI
  fprintf(['writeCTFMRI: Structure array MRItags does not contain the tags ',...
      'for a CTF MRI file (v4 format).\n'],fileName);
  return
end
clear isCTFMRI

%  Add clinical use and creator software messgaes
MRItags=update_MRI_descriptors(MRItags,clinical_use_message,creatorSoftware);

%  Add MRIdata to MRItags, and add start and end tags
nTag=length(MRItags);
MRItags(2:nTag+1)=MRItags(1:nTag);
MRItags(1)=struct('name','WS1_','type',0,'data',[]);  % Start tag
nTag=nTag+1;
for k=1:nSlice
  MRItags(nTag+k)=struct('name',['_CTFMRI_SLICE_DATA#',num2str(k,'%5.5d')],'type',3,...
    'data',int16(reshape(MRIdata(:,:,k),MRIsize^2,1)));
end
nTag=nTag+nSlice;
MRItags(nTag+1)=struct('name','EndOfParameters','type',-1,'data',[]);  % End tag
clear MRIdata k nSlice MRIsize nTag;

if isempty(printWarning)
  fprintf(['\nwriteCTFMRI: The MRI you are creating has been processed by software not\n',...
      '\tmanufactured by VSM MedTech Ltd. and that has not received marketing clearance\n',...
      '\tfor clinical applications.  This MRI should not be later employed for clinical\n',...
      '\tand/or diagnostic purposes.\n']);
  printWarning=1;
end

%  Create the MRI file.
writeCPersist(fileName,MRItags);
return
%%%%%%%%%%%%%%  end of readMRI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Function update_MRI_descriptors   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MRItags=update_MRI_descriptors(MRItags,comment,creatorSoftware);

%  Makes sure that certain tags are in structure MRItags.

%  Inputs: MRItags : Structure array produced by readCTFMRI.
%          comment : Character string that is added to infods tags listed in addCommentTag.
%          creatorSoftware : Character string indicating that the MRI set was
%              created by writeCTFMRI.  Added to the tags listed in addCreatorTag.

%  Adds comment (clinical use message) and creator software name to 

%  MRI tags that will display the comment.
addCommentTag=strvcat('_PATIENT_ID','_PATIENT_NAME',...
  '_STUDY_COMMENTS','_CTFMRI_COMMENT');

%  MRI tags that will display the creator software.
addCreatorTag=strvcat('_CTFMRI_COMMENT');

%  Check that the comment and creator tags are in structure MRItags.  If not, issue an
%  error message and set MRItags to [].
addTag=strvcat(addCommentTag,addCreatorTag);
addIndex=zeros(1,size(addTag,1));
k=0;
for k=1:length(MRItags)
  tagIndex=strmatch(MRItags(k).name,addTag,'exact');
  addIndex(tagIndex)=k;
  if all(addIndex)~=0;break;end
end
%  addTag(j) is MRItags(addIndex(j))

if any(addIndex==0)
  %  At least one of the addTags is missing.  Print an error message and return.
  fprintf(['writeCTFMRI (update_MRI_descriptors): At least one tag is missing from ',...
      'structure array MRItags.\n']);
  for k=1:size(addTag,1)
    if isempty(strmatch(addTag(k,:),addTag(1:k-1,:),'exact'))
      fprintf('\t\t\t%s\n',deblank(addTag(k,:)));
    end
  end
  MRItags=struct([]);  % Force an error in the calling program.
  return
end

%  Check that all of the addTags are text strings
type=[MRItags.type];
if any(type(addIndex)~=9 & type(addIndex)~=10)
  fprintf('writeCTFMRI (update_MRI_descriptors): At least one of the new tags is not a character string.\n');
  for k=1:length(addTag)
    if isempty(strmatch(addTag(k,:),addTag(1:k-1,:),'exact'))
      fprintf('\t\t\t%s   type=%d\n',deblank(addTag(k,:)),type(addIndex(k)));
    end
  end
  MRItags=struct([]);  % Force an error in the calling program.
  return
end

addCommentIndex=addIndex(1:size(addCommentTag,1));
addCreatorIndex=addIndex(size(addCommentTag,1)+[1:size(addCreatorTag,1)]);

if exist('creatorSoftware')~=1;
  creatorSoftware=char([]);
elseif ~ischar(creatorSoftware)
  creatorSoftware=char([]);
else
  creatorSoftware=deblank(creatorSoftware);
end

if exist('comment')~=1;
  comment=char([]);
elseif ~ischar(comment)
  comment=char([]);
else
  comment=deblank(comment);
end

for k=addCreatorIndex
  if isempty(MRItags(k).data)
    MRItags(k).data=creatorSoftware;
  else
    MRItags(k).data=[creatorSoftware '  ' MRItags(k).data];
  end
end
  
for k=addCommentIndex
  if isempty(MRItags(k).data)
    MRItags(k).data=comment;
  else
    MRItags(k).data=[comment '  ' MRItags(k).data];
  end
end

%  Shorten type 9 (CStr32) strings to 31 characters plus a terminating null.
for k=addIndex
  if MRItags(k).type==9 & length(MRItags(k).data)>31;
    MRItags(k).data=[MRItags(k).data(1:31) char(0)];
  end
end

return
%%%%%%%%%%%%%% End of  update_MRI_descriptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
