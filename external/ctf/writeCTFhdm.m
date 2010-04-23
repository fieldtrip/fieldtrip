function writeCTFhdm(fileName,hdm);

%  Version 1.1   20 April 2007   Modified to make sure that clinical-0ise message and 
%                                creator software are not repeated.
%                                Test date.
%                22 March 2007.  Modified to write v6.0 .hdm files which have additional
%                                fields in MultiSphere_Data.  Check that the .hdm version
%                                number is compatible with the fields included in
%                                MultiSphere_Data.
%  Creates a CTF-compatible head model file from structure hdm
%  Structure hdm is in the format produced by readCTFhdm.

%  The formal of hdm (Head Model) files is given in document "CTF MEG FIle Formats",
%  PN900-0088.  The head model files are text files in the "Config Reader" format.

%  The purpose is to allow MATLAB users to specify multiple-local spheres head models
%  in a formatthat is compatible with CTF analysis software (DipoleFit).  It also allows
%  users to transfer .hdm files between datasets.

%  Head Model File format is defined in document "CTF MEG FIle Formats', PN900-0088.

%  Inputs : fileName: Name of the .hdm file to create.  It must include the full path,
%                     including dataset specification and the extension .hdm.
%           hdm:  The structure containing the class information.  See D

%  Output: Head model file.

persistent printWarning
creatorSoftware='writeCTFhdm';
versionSoftware='v1.1';

if nargin==0
  fprintf('writeCTFhdm   Version 1.1   20 April 2007\n');
  fprintf(['\twriteCTFhdm(fileName,hdm) creates a CTF-format head model file',...
      ' from the fields of structure hdm.\n',...
      '\t\t\tfileName is the name of the new head model file including path and ',...
      'extension ''.hdm''.\n',...
      '\t\t\thdm is a structure with the format produced by readCTFhdm.\n',...
      '\t\t\tThe head model file format is described in document ',...
      '"CTF MEG File Formats", PN900-0088.\n\n']);
  return
elseif nargin~=2
  fprintf(['\nwriteCTFhdm:  %d input argument(s).  writeCTFhdm needs 2 arguments:\n',...
      '              fileName and structure containing the head model.\n\n'],nargin);
  return
elseif ~ischar(fileName) | isempty(fileName) | ~isstruct(hdm) | isempty(hdm)
  fprintf('writeCTFhdm:  Wrong argument types, or empty arguments.\n');
  whos fileName hdm;
  return
elseif ndims(fileName)>2 | min(size(fileName))>1
  fprintf('\nwriteCTFhdm: size(fileName)=[');fprintf(' %d',size(fileName));
  fprintf(']    Input fileName must be a character string.\n\n');
  return
elseif exist(fileName)==2
  fprintf('writeCTFhdm: File %s already exists.\n',fileName);
  return
elseif isempty(strfind(fileName,'.hdm'))
  fprintf('writeCTFhdm: fileName=%s does not have .hdm extension.\n',fileName);
  return
elseif ~isfield(hdm,'Model') | ...
    (~isfield(hdm,'MultiSphere_Data') & ~isfield(hdm,'MEG_Sphere'))
  fprintf(['writeCTFhdm: Structure hdm is missing fields Multisphere_Data, ',...
      'MEG_Sphere or Model.\n']);
  fprintf('             It cannot be a useful head model file.\n');
  hdm
  return
elseif isfield(hdm,'MEG_Sphere') & ~isfield(hdm,'MultiSphere_Data')
  fprintf(['writeCTFhdm: Creating a single-sphere model.  ',...
      'No multisphere data found in structure hdm.\n']);
elseif ~isfield(hdm,'MEG_Sphere') & isfield(hdm,'MultiSphere_Data')
  fprintf(['writeCTFhdm: Creating a multi-sphere model.  ',...
      'No single-sphere data found in structure hdm.\n']);
end

%  Check version numbers and 'MultiSphere_Data fields.
if isfield(hdm,'File_Info') & isfield(hdm,'MultiSphere_Data')
  versionText=hdm.File_Info.VERSION;
  kpt=strfind(versionText,'VERSION_');
  version=sscanf(versionText(kpt+8:length(versionText)),'%g');
  if version<=5.9
    if isfield(hdm.MultiSphere_Data,'HEADPOS') | isfield(hdm.MultiSphere_Data,'SURFACE_TYPE')
      fprintf('\n');
      fprintf('%s\n',['writeCTFhdm: .hdm VERSION=',num2str(version,'%0.1f'),...
          ', but fields HEADPOS and/or SURFACE_TYPE appear'],...
        '              in hdm.MultiSphere_data.  They should appear only in v6.0 and higher.');
      return
    end
  elseif version>=6.0
    if ~isfield(hdm.MultiSphere_Data,'HEADPOS') | ~isfield(hdm.MultiSphere_Data,'SURFACE_TYPE')
      fprintf('\n');
      fprintf('%s\n',['writeCTFhdm: .hdm  VERSION=',num2str(version,'%0.1f'),...
          ', but fields HEADPOS and/or SURFACE_TYPE do not appear'],...
        '              in hdm.MultiSphere_data.  They must appear in v6.0 and higher.');
      return
    end
  end
end

if isempty(printWarning)
  fprintf(['\nwriteCTFhdm: The head-model data you are writing have been processed by software\n',...
      '\tnot manufactured by VSM MedTech Ltd. and that has not received marketing clearance\n',...
      '\tfor clinical applications.  These data should not be later employed for clinical\n',...
      '\tand/or diagnostic purposes.\n']);
  printWarning=1;
end

fid=fopen(fileName,'w','ieee-be');  % Use 'w' option to ensure no char(13)'s at the end of line.

%  Add the clinical-use message if it is not already present.
addClinicalUseMessage=1;
addCreatorSoftwareMessage=1;
if isfield(hdm,'note')
  addClinicalUseMessage=isempty(strfind(reshape(hdm.note',1,prod(size(hdm.note))),'NOT FOR CLINICAL USE'));
  addCreatorSoftwareMessage=isempty(strfind(reshape(hdm.note',1,prod(size(hdm.note))),'writeCTFhdm'));
end
if isfield(hdm,'note')
  for k=1:size(hdm.note,1)
    strng=deblank(hdm.note(k,:));
    if ~strcmp(strng(1:2),'//');strng=['// ',strng];end
    fprintf(fid,'%s\n',strng);
    if k==3 & addClinicalUseMessage;
      fprintf(fid,[...
          '// *************************************\n',...
          '// *    NOT FOR CLINICAL USE           *\n',...
          '// *************************************\n\n']);
    end
  end
  clear strng;
end
if addCreatorSoftwareMessage
  fprintf(fid,'//\n// Prepared by %s %s\n//\n',creatorSoftware,versionSoftware);
end
if addClinicalUseMessage
  fprintf(fid,['// *************************************\n',...
      '// *    NOT FOR CLINICAL USE           *\n',...
      '// *************************************\n\n']);
end
clear addClinicalUseMessage;

className=char(fieldnames(hdm));
for n=1:size(className,1)
  cName=deblank(className(n,:));
  if strcmp(cName,'note');
    continue;  % Already printed the notes.
  end
  eval(['tagName=char(fieldnames(hdm.',cName,'));']);
  fprintf(fid,'%s\n{\n',cName);
  if ~strcmp('MultiSphere_Data',cName);
    for k=1:size(tagName,1)
      tgName=deblank(tagName(k,:));
      if strcmp(tgName,'note');
        eval(['note=hdm.',cName,'.note;']);
        if ~isempty(note)
          for q=1:size(note,1)
            strng=deblank(note(k,:));
            if ~strcmp(strng(1:2),'//');strng=['// ',strng];end
            fprintf(fid,'\t%s\n',strng);
          end
        end
        clear strng;
      else
        eval(['tgVal=hdm.',cName,'.',tgName,';']);
        tgType=class(tgVal);
        fprintf(fid,'\t%s:',tgName);
        if isempty(tgVal);
          fprintf(fid,'\t');
        elseif strcmp(tgType,'char')
          fprintf(fid,'\t%s',tgVal);
        elseif max(abs(tgVal-round(tgVal)))<1e-4
          fprintf(fid,'\t%d',round(tgVal));
        else
          fprintf(fid,'\t%0.4f',tgVal);
        end
        fprintf(fid,'\n');
        clear tgType tgVal tgName;
      end
    end
  else  % MultiSphere_Data has a different format :
    if strmatch('SEARCH_RADIUS',tagName,'exact')
      fprintf(fid,'\tSEARCH_RADIUS:\t%5.3f\n',hdm.MultiSphere_Data.SEARCH_RADIUS);
    else
      fprintf('writeCTFhdm: tag SEARCH_RADIUS is missing from class MultiSphere_Data.\n');
    end
    if strmatch('HEADSHAPE_FILE',tagName,'exact')
      fprintf(fid,'\tHEADSHAPE_FILE:\t%s\n',hdm.MultiSphere_Data.HEADSHAPE_FILE);
    else
      fprintf('writeCTFhdm: tag HEADSHAPE_FILE is missing from class MultiSphere_Data.\n');
    end
    if strmatch('SURFACE_TYPE',tagName,'exact')
      fprintf(fid,'\tSURFACE_TYPE:\t%s\n',hdm.MultiSphere_Data.SURFACE_TYPE);
    end
    fprintf(fid,'\n');
    if strmatch('HEADPOS',tagName,'exact')
      fprintf(fid,'\t%s\n','// Head Coil coordinates relative to dewar, cm',...
        '// 			NA_x	NA_y	NA_z	LE_x	LE_y	LE_z	RE_x	RE_y	RE_z	Nominals');
      fprintf(fid,'\tHEADPOS:');
      fprintf(fid,'\t%0.3f',hdm.MultiSphere_Data.HEADPOS.HEADPOS);
      fprintf(fid,'\t%s\n\n',hdm.MultiSphere_Data.HEADPOS.NOMINAL);
    end
    
    %  Don't use hdm.MultiSphere_data.note
    fprintf(fid,'\t%s\n','// Multiple Sphere locations in cm');
    fprintf(fid,'\t%s','//','X','Y','Z','Radius');
    fprintf(fid,'\n');
    %  Create a separate tag for each SQUID channel.  Be carefule to remove the 
    %  coefficient-set identifier in case the user has generated array SQUIDname from
    %  ds.res4.chanNames of some dataset.
    for q=1:size(hdm.MultiSphere_Data.SQUIDname,1)
      fprintf(fid,'\t%s:\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n',...
        deblank(strtok(hdm.MultiSphere_Data.SQUIDname(q,:),'-')),...
        hdm.MultiSphere_Data.sphereOrigin(:,q),...
        hdm.MultiSphere_Data.radius(q));
    end
    fprintf(fid,'\n');
  end
  fprintf(fid,'}\n\n');
  clear cName tagName k q;
end  % End of loop over classes
fclose(fid);
return