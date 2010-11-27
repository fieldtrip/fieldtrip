%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Function writeCPersist  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeCPersist(filename,Tag)

%  Version 1.2 24 April 2007.  Removed terminating null characters from charatcer strings
%                              types 10,11)

%  Version 1.1   13 April 2007
%  Write a CTF CPersist file.  See document CTF MEG File Format, PN900-0066  

%  Inputs : filename : Name of the output file name including path
%           Tag: Structure array containing all of the tags, types, and data.

% 31 Aug 2006: Recognizes type=1,...,17.  If other data types are specified,
%              writeCPersist will print an error message.

%  Delete existing file so the new file has the correct creation data/time.

if nargin==0
  fprintf(['writeCPersist: Version 1.2   24 April 2007    Writes a CTF CPersist file.\n'...
      '\twriteCPersist(filename,Tag) writes a CPersist file from the contents of\n',...
      '\tstructure array Tag.  Tag is in the format prepared by readCPersist.\n\n']);
  Tag=[];
  return
end

if exist(filename)==2
  delete(filename);
end

startString='WS1_';
EOFstring='EndOfParameters';

startVal=256.^[3:-1:0]*double(startString)';  %  Integer version of startString

fid=fopen(filename,'w','ieee-be');

EOFcount=0;  % Make sure that the startString's and EOFstring's balance
count=0;
while count<length(Tag)
  count=count+1;
  if strcmp(Tag(count).name,startString)     % start of Cpersist object
    fwrite(fid,startVal,'int32');
    EOFcount=EOFcount-1;
    continue
  end
  fwrite(fid,length(Tag(count).name),'int32'); % end of CPersist object
  fwrite(fid,Tag(count).name,'char');
  if strcmp(Tag(count).name,EOFstring);
    EOFcount=EOFcount+1;
    continue;
  end
  fwrite(fid,Tag(count).type,'int32');
  if Tag(count).type==1
    fwrite(fid,Tag(count).data,'int32');
  elseif Tag(count).type==2   % Start embedded CPersist object
  elseif Tag(count).type==3   % binary list 
    fwrite(fid,2*length(Tag(count).data),'int32');
    fwrite(fid,Tag(count).data,'int16');
  elseif Tag(count).type==4   % double
    fwrite(fid,Tag(count).data,'float64');
  elseif Tag(count).type==5    % integer
    fwrite(fid,Tag(count).data,'int32');
  elseif Tag(count).type==6  % short integer
    fwrite(fid,Tag(count).data,'int16');
  elseif Tag(count).type==7  % unsigned short integer
    fwrite(fid,Tag(count).data,'uint16');
  elseif Tag(count).type==8  % Boolean byte
    fwrite(fid,Tag(count).data,'uint8');
  elseif Tag(count).type==9;   % CStr32
    nChar=min(32,length(Tag(count).data));
    fwrite(fid,[double(Tag(count).data(1:nChar)) zeros(1,32-nChar)],'uint8');
  elseif Tag(count).type==10;   % CString
    fwrite(fid,length(Tag(count).data),'int32');
    fwrite(fid,double(Tag(count).data),'uint8');
  elseif Tag(count).type==11;  % Cstring list.
    nList=size(Tag(count).data,1);
    fwrite(fid,nList,'int32');
    for k=1:nList
      %  Do not force termination of strings with nulls (char(0))
      Cstring=[deblank(Tag(count).data(k,:))];
      fwrite(fid,length(Cstring),'int32');
      fwrite(fid,double(Cstring),'uint8');
    end
    clear k CString nList;
  elseif Tag(count).type==12;  % CStr32 list.
    nList=size(Tag(count).data,1);  %size(data)=[nList 32]
    fwrite(fid,nList,'int32');
    %  Do not force termination of strings with nulls (char(0))
    for k=1:nList
      strng=deblank(Tag(count).data(k,:));
      nChar=min(32,length(strng));
      Tag(count).data(k,:)=[strng(1:nChar) char(zeros(1,32-nChar))];
    end
    fwrite(fid,double(Tag(count).data)','uint8');
    clear k strng nChar nList;
  elseif Tag(count).type==13;  % SensorClass list.
    fwrite(fid,length(Tag(count).data),'int32');
    fwrite(fid,Tag(count).data,'int32')
  elseif Tag(count).type==14    % long integer
    fwrite(fid,Tag(count).data,'int32');
  elseif Tag(count).type==15    % unsigned longinteger
    fwrite(fid,Tag(count).data,'uint32');
  elseif Tag(count).type==16    % unsigned integer
    fwrite(fid,Tag(count).data,'uint32');
  elseif Tag(count).type==17    % Boolean
    if ~any(Tag(count).data==[0 1])
      fprintf('writeCPersist: tagname=%s type=%d value=%d?  (Must =0 or 1)\n',...
        tagname,type,Tag(count).data);
      Tag(count).data=(Tag(count).data~=0);
      fprintf('               Set value=%d.\n',Tag(count).data);
    end
    fwrite(fid,Tag(count).data,'int32');
  else
    fprintf('writeCPersist: Tag(%d) name=%s  type=%d\n',count,Tag(count).name,Tag(count).type);
    fprintf('\t\t UNRECOGNIZED type.\n');
    fclose(fid);
    break;
  end
end   % Loop over tags
fclose(fid);

% EOF_count should be zero at the end of the file
if EOFcount~=0;
  fprintf('write_CPerist: EOFcount=%d   Stop strings do not balance start strings.\n',...
    EOFcount);
  fprintf('               Start string=%s    Stop string=%s\n',startString,EOFstring);
end
return
%%%%%%%%  End of writeCPersist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
