function [ret]=ml_GetPrivateProfileString(AppName, KeyName, filename)
% Usage: [ret]=ml_GetPrivateProfileString(AppName, KeyName, filename)
% This function is the matlab simplification of the Private Profile String
% Reading utilities in VB and C++ to handle .ini file data
% AppName = Section name of the type [xxxx]
% KeyName = Key name separated by an equal to sign, of the type abc = 123
% filename = ini file name to be used
% ret = string value of the key found, 
% 
% Example: a sample ini file (sample.ini) may have entries:
%       [XYZ]
%       abc = 123
%       [ZZZ]
%       abc = 890
% ml_GetPrivateProfileString('XYZ','abc', 'sample.ini') will return '123' 
% ml_GetPrivateProfileString('ZZZ','abc', 'sample.ini') will return '890' 
% if key or the appname is not found, it returns '** Error **' string
% if the key is empty, then an empty string is returned
%
% Created By: Irtaza Barlas
% Created On: June 9, 2004
% Created For: IAS Inc.

try
    ret = '** Error **';
a=0;
if nargin~=3
    return
end
if ~exist(filename,'file')
    return
end
if isempty(AppName)==1 || isempty(KeyName)==1 
    return
end
a = fopen(filename);
if a<=0
    return
end

appstr=['[' AppName ']'];
found_app=0;
st = fgetl(a);
while found_app~=1
    if isempty(st)==1 
        st = fgetl(a);
        continue
    elseif st==-1
        break
    end
   if strcmp(st, appstr)>0
       % look for the key
       found_app=1;
       found_key=0;
       kst = fgetl(a);
       while found_key==0
           if isempty(kst)==1
               kst = fgetl(a);
               continue
           end
           if kst==-1
               break
           end
           if isempty(kst) == 0 
               kst = sscanf(kst, '%s'); % helps eliminate whitespaces
               if kst(1)=='['   % next key                   
                   break
               end
               % find equal to sign
               eq_idx=find(kst=='=');
               if isempty(eq_idx)~=1 && eq_idx(1)>1 
                   key=kst(1:eq_idx(1)-1);
                   if strcmp(key, KeyName)>0 
                       [rs, cs]=size(kst);
                       if eq_idx(1)>=cs
                           ret = '';
                       else
                           ret=kst(eq_idx(1)+1:end);
                       end
                       found_key=1;
                       break
                   end
               end               
           end
           kst = fgetl(a);
       end
       break
   end
   st = fgetl(a);
end
fclose(a);
catch
    if a>0 
        fclose(a);
    end
    [ret, rid] = lasterr;
end
