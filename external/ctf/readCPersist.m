function Tag=readCPersist(CPfile,printtags)

%  Version 1.2   24 April 2007   Modified to close the CPersist file if a really huge
%  taglength is encountered.  Recently discovered .acq files with the string 'ssss' added
%  at the end of the file after the final 'EndofParameters' string.

%  Reads a CTF CPersist file.  See document CTF MEG File Format, PN900-0066  
%  printtags=1 : Print a list of tags and values.  If not supplied, default=0;
%  Recognizes data types 1,...,17.  If other types appear in the file, an error
%  message is printed.

%  Input : CPfile: Name of a CPersist file including the full path.
%          printtags: 0 : Do not print the tags as they are read.
%                     1 : Do print the tags.

%  Output : Tag : Structure array containing all of the tags, types, and data.
%
%  Adds 2 tag types: 0 marks start string('WS1_')
%                   -1 marks 'EndofParameters'

EOPcountMessage=0;
%  A reading error could give a huge tagLength which would result in an
%  out-of-memory error when executing command tagName=fread(fid,tagLength,'uint8');
%  Limit acceptable values of tagLength so this will not happen.
maxNameLength=2000;  % Do not accept tags with names exceeding maxNameLength characters.

if nargin==0
  fprintf(['readCPersist: Version 1.2   24 April 2007    Reads a CTF CPersist file.\n\n'...
      '\tTag=readCPersist(filename) or Tag=readCPersist(filename,0) reads\n',...
      '\ta CPersist file and stores the contents in structure array Tag.\n\n',...
      '\tTag=readCPersist(filename,1) reads a CPersist file and prints ',...
      'the contents.\n\n',...
      '\tTag types are defined in document CTF MEG File Formats, PN900-0088,',...
      ' Appendix B.\n',...
      '\treadCPersist also creates type=0 tags to mark the beginning and type=-1\n',...
      '\ttags to mark the end of CPersist objects within the file.\n\n']);
  if nargout>1;Tag=[];end
  return
end

startString='WS1_';
EOPstring='EndOfParameters';

startVal=256.^[3:-1:0]*double(startString)';

% try with ending delimiter to dereference symlinks
D=dir([CPfile,filesep]);
% if empty, try accessing the file directly
if isempty(D)
  D=dir(CPfile);
end

if isempty(D)
  fprintf('readCPersist:  Could not find file %s.\n',CPfile);
  Tag=struct([]);
  return
elseif D.bytes==0
  fprintf('readCPersist:  File %s contains 0 bytes.\n',CPfile);
  Tag=struct([]);
  return
end

totalBytes=D.bytes;
clear D;

if ~exist('printtags');printtags=[];end
if ~isnumeric(printtags) | ~isequal(size(printtags),[1 1]);printtags=0;end

fid=fopen(CPfile,'r','ieee-be');

%  Step through tags and info.
tagCount=0;
bytePtr=0;  % No. of bytes read so far.  
% Compare bytePtr to totalBytes to decide if an impossible read request is being made.
% This is necessary because some .acq files have junk string 'ssss' appended after the
% final EOPstring.
EOPcount=0;  % Counts end-of-object markers.  
%  EOPcount is incremented at the start of an object (when startString is read).
%  EOP is decremented at the end of an object (EOPstring is read).

while 1
  tagCount=tagCount+1;
  %  tagLength = number of characters in the tagName.
  tagLength=fread(fid,1,'int32');
  bytePtr=bytePtr+4;
  if isempty(tagLength) | bytePtr>=(totalBytes-8)
    break;
  elseif tagCount==1 & tagLength~=startVal  % Is this a CPersist file?
    fprintf(['readCPersist: File %s does not start with %s',...
        ' (required for a CPersist file).\n'],CPfile,startString);
    fclose(fid);
    Tag=[];
    return
  elseif tagLength==startVal   % Start of embedded CPersist object
    tagName=startString;
    buff=[];
    type=0;
    EOPcount=EOPcount+1;
    if printtags;
      fprintf('%4d: %s  EOPcount=%d  bytePtr=%5d',tagCount,tagName,EOPcount,bytePtr);
    end
  elseif tagLength<=0 | tagLength>maxNameLength | (bytePtr+tagLength)>totalBytes
    %  Should never come here, but some versions of Acq add junk bytes after the final
    %  EOPstring, and they may be interpreted as the number of  bytes in the next tag.
    if tagLength<0
      fprintf('readCPersist: tagCount=%d  bytePtr=%d  tagLength=%d??\n',tagCount,bytePtr,tagLength);
    end
    break;
  else
    %  Read the tagName
    tagName=char(fread(fid,tagLength,'uint8'))';
    bytePtr=bytePtr+tagLength;
    if length(tagName)~=tagLength
      fprintf(['readCPersist: tagCount=%d  tagLength=%d  tagName=%s',...
          ' (length~=tagLength)\n'],tagCount,tagLength,tagName);
      break;
    elseif strcmp(EOPstring,tagName)   % end of CPersist object
      EOPcount=EOPcount-1;
      if EOPcountMessage==0 & EOPcount<0 & printtags==0
        %  Error in the file structure, but keep reading because it may help user to see
        %  partial results from reading a bad file.
        fprintf('%s\n',' ',['readCPersist:  More end-of-object marks (',EOPstring,...
            ') than start-of-object marks (',startString,') have been read so far.'],...
          ['               File ',CPfile,' has a format error.'],' ');
        EOPcountMessage=1;
      end
      buff=[];
      type=-1;
      if printtags~=0;
        fprintf('%4d:  %s                 EOPcount=%d',tagCount,tagName,EOPcount);
      end
    else
      type=fread(fid,1,'int32');
      bytePtr=bytePtr+4;
      if isempty(type)
        fprintf('readCPersist: tagCount=%d tagName=%s  type=[]??\n',tagCount,tagName);
        break;
      end
      if printtags;fprintf('%4d:  %s  type=%d',tagCount,tagName,type);end
      if type==1      % custom type  Assume zero length except for DatasetFiles
        if strcmp(tagName,'DatasetFiles');
          buff=fread(fid,1,'int32');
          bytePtr=bytePtr+4;
        else
          buff=[];
        end  
      elseif type==2  % nested CPersist object
        buff=[];
      elseif type==3  % binary list 
        bufflength=fread(fid,1,'int32');  % buffer in bytes.
        bytePtr=bytePtr+4;
        if bufflength~=2*round(bufflength/2)
          fprintf(['readCPersist:  ERROR : tagName=%s  bufflength=%d',...
              '  MUST BE EVEN.\n'],round(bufflength));
          break;
        end    
        buff=int16(fread(fid,round(bufflength/2),'int16')');
        bytePtr=bytePtr+bufflength;
        if printtags;fprintf('  %d entries',round(bufflength/2));end
      elseif type==4    % double
        buff=fread(fid,1,'float64');
        bytePtr=bytePtr+8;
        if printtags;fprintf('  %f',buff);end
      elseif type==5    % integer
        buff=fread(fid,1,'int32');
        bytePtr=bytePtr+4;
        if printtags;fprintf('  %d',buff);end
      elseif type==6    % short integer
        buff=fread(fid,1,'int16');
        bytePtr=bytePtr+2;
        if printtags;fprintf('  %d',buff);end
      elseif type==7    % unsigned short integer
        buff=fread(fid,1,'uint16');
        bytePtr=bytePtr+2;
        if printtags;fprintf('  %d',buff);end
      elseif type==8    % Boolean byte
        buff=char(fread(fid,1,'uint8')');
        bytePtr=bytePtr+1;
        if printtags;fprintf('  %d',buff);end
      elseif type==9    % CStr32  32-character string
        fprintf('\ttype=9  bytePtr=%d  tagLength=%d\n',bytePtr,tagLength);
        buff=char(fread(fid,32,'uint8')');
        bytePtr=bytePtr+32;
        if printtags;fprintf('   %s',buff);end
      elseif type==10;   % CString
        %fprintf('\ttype=10  bytePtr=%d  tagLength=%d\n',bytePtr,tagLength);
        stringlength=fread(fid,1,'int32');
        bytePtr=bytePtr+4;
        buff=char(fread(fid,stringlength,'uint8')');
        bytePtr=bytePtr+stringlength;
        if printtags;
          fprintf('  %d chars',stringlength);
          if stringlength<80;fprintf('   %s',buff);end
        end
      elseif type==11;  % Cstring list.
        nList=fread(fid,1,'int32');
        %fprintf('\ttype=11  bytePtr=%d  tagLength=%d  nList=%d\n',bytePtr,tagLength,nList);
        bytePtr=bytePtr+4;
        buff=[];
        for k=1:nList
          nLen=fread(fid,1,'int32');
          buff=strvcat(buff,char(fread(fid,nLen,'uint8')'));
          bytePtr=bytePtr+4+nLen;
        end
        if printtags;fprintf('  List of %d strings.',nList);end
      elseif type==12;  % CStr32 list.
        nList=fread(fid,1,'int32');
        fprintf('\ttype=12  bytePtr=%d  tagLength=%d  nList=%d\n',bytePtr,tagLength,nList);
        buff=reshape(char(fread(fid,32*nList,'uint8')),32,nList)';  %size(buff)=[nList 32]
        bytePtr=bytePtr+4+32*nList;
        if printtags;fprintf('  List of %d strings.',nList);end
      elseif type==13;  % SensorClass list.
        nList=fread(fid,1,'int32');
        buff=fread(fid,nList,'int32');
        bytePtr=bytePtr+4+4*nList;
        if printtags;
          fprintf('  %d-element list.',nList);
          if nList<10;
            fprintf(' %d',buff);
          end
          fprintf('\n');
        end
      elseif type==14    % long integer
        buff=fread(fid,1,'int32');
        bytePtr=bytePtr+4;
        if printtags;fprintf('  %d',buff);end
      elseif type==15    % unsigned longinteger
        buff=fread(fid,1,'uint32');
        bytePtr=bytePtr+4;
        if printtags;fprintf('  %d',buff);end
      elseif type==16    % unsigned integer
        buff=fread(fid,1,'uint32');
        bytePtr=bytePtr+4;
        if printtags;fprintf('  %d',buff);end
      elseif type==17    % Boolean
        buff=fread(fid,1,'int32');
        bytePtr=bytePtr+4;
        if buff~=0 & buff~=1
          fprintf('readCPersist: tagName=%s type=%d value=%d?  (Must =0 or 1)\n',...
            tagName,type,buff);
          buff=(buff~=0);
          fprintf('               Set value=%d.\n',buff);
        end
      else
        %  Should never come here.
        fprintf('readCPersist:  tagName=%s type=%d',tagName,type);
        fprintf('\t\t UNRECOGNIZED type.\n');
        break;
      end
    end
  end
  %  Assemble the next element in the structure array.
  Tag(tagCount)=struct('name',tagName,'type',type,'data',buff);
  clear tagName type buff;
  if printtags;fprintf('\n');end
end
fclose(fid);

% EOPcount should be zero at the end of the file
if EOPcount~=0;
  fprintf(['readCPersist: Stop strings (',EOPstring,') do not balance start strings (',...
      startString,').\n              (No. of start strings)-(No. of stop strings)=',...
      int2str(EOPcount),'.\n']);
end
return
