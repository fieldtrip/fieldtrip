function [err] = TellAfni (cs, opt)
%
%   [err] = TellAfni (cs, opt)
%
%Purpose:
%   Drive AFNI
%   
%   
%Input Parameters:
%   cs: An Nx1 vector of communication command structures
%       It is obtained using NewCs structure
%   opt: An optional options structure
%        .QuitOnErr (0/[1]): Return from function if any of cs(i) is malformed
%        .Verbose (0/[1]/2): 0 = mute, 1 = Yak, Yak (default), 2 = YAK YAK
%Output Parameters:
%   err : 0 No Problem
%       : N  N Problems
%   
%      
%More Info :
%   
%     NewCs
%     TellAfni_Commands
%     Test_TellAfni
%     AFNI's README.driver  file and the program plugout_drive
%     
%     Author : Ziad Saad
%     Date : Tue Dec 6 10:38:23 EST 2005
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'TellAfni';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

LogFile = '.TellAfni.log';

if (nargin == 1) opt.QuitOnErr = []; end

if (~isfield(opt, 'QuitOnErr') | isempty(opt.QuitOnErr)) opt.QuitOnErr = 1; end
if (~isfield(opt, 'Verbose') | isempty(opt.Verbose)) opt.Verbose = 1; end

ncs = length(cs);
if (ncs == 0) err = 0; return; end

if (isempty(LogFile) | LogFile(1) ~= '.'),
   fprintf(2,'Bad LogFile, whas happening?');
   return;
end

com = '';
ComLast = '';
ncom = 0;
for (i=1:1:ncs),
   if (cs(i).err),
      fprintf(2,'Error in command structure %d (%s).\n', i,cs(i).c);
      if (opt.QuitOnErr),
         fprintf(2,'Quitting.\n');
         return;
      else
         fprintf(2,'Ignoring command.\n');
         cs(i).c = '';
      end
   end
   
   if (cs(i).c),
      switch (cs(i).c),
         case 'START_AFNI',
            if (exist(LogFile) == 2),
               fprintf(2,'\nError: Log file %s found.\n', LogFile);
               fprintf(2,'Make sure no other AFNI is running AND listening for plugouts.\n');
               fprintf(2,'rm the logfile with:\n!rm -f %s \n', LogFile);
               fprintf(2,'Then try your command again.\n\n');
               return;
            end
            scom = sprintf('afni -yesplugouts %s |& tee %s &', cs(i).v, LogFile);
            [s,w] = unix(scom);
            if (opt.Verbose & ~isempty(w)), 
               fprintf(1,'Command output:\n%s\n', w);
            end
            if (s),
               fprintf(2,'Error launching afni\n');
               return;
            end
            %check on errors regarding connections
            iserr = 1; Tel = 0;
            while (Tel < 10 & iserr),
               if ( (rem(Tel,1) < 0.0001) ) fprintf(2,'Checking on communication (%.1f/10)\n', Tel); end
               iserr = TellAfni_CheckLog(LogFile, (rem(Tel,1)<0.0001));
               pause (0.2); % wait to be sure AFNI launched
               Tel = Tel+0.2;
            end
            if (iserr),
               fprintf(2,'\nLaunched a new AFNI but failed to communicate.\n');
               fprintf(2,'Close the failed AFNI and any other AFNI sessions\n');
               fprintf(2,'that are listening to plugouts.\n\n');
               scom = sprintf('rm -f %s >& /dev/null',  LogFile);
               unix(scom);
               return;
            end
         case 'QUIT',
            if (~isempty(LogFile)),
               ComLast = sprintf('rm -f %s >& /dev/null', LogFile);
            end
            com = sprintf('%s -com ''%s %s''', com, cs(i).c, cs(i).v);
            ncom = ncom+1;
         otherwise,
            com = sprintf('%s -com ''%s %s''', com, cs(i).c, cs(i).v);
            ncom = ncom+1;
      end
   end
end

b = 0;
if (~isempty(com)),
   scom = sprintf('plugout_drive -v %s -quit', com);
   if (opt.Verbose > 1),
      fprintf(1, 'making call:\n%s\n', scom);
   end
   [s,w] = unix(scom);
   if (isempty(strfind(scom,'QUIT'))), %Do not check in cases of QUIT
      [err, g, b] = TellAfniCheck(w);
      if (ncom ~= err + g + b ),
         fprintf(2,'Warning: Unexpected parsing trouble (ncom=%d, err+g+b=%d).\n', ncom, err+g+b);
      end
      if (err),
         fprintf(2,'Warning: Failed in parsing plugout_drive output.\nCannot confirm how %d out of %d commands executed\n', err, ncom);
      end
      if (g & opt.Verbose),
         fprintf(1,'%d out of %d commands OK\n', g, ncom);
      end
      if (b),
         fprintf(2,'Warning: %d out of %d commands failed in AFNI.\n', b, ncom);
      end
   end
   if (opt.Verbose > 1), 
      fprintf(1,'Command output:\n%s\n', w);
   end   
   if (s),
      fprintf(2,'Error telling afni\n');
      return;
   end
end

if (~isempty(ComLast)),
   [s,w] = unix(ComLast);
end

err = b;
return;

function err = TellAfni_CheckLog(LogFile, verb)

err = 1;
   if (exist(LogFile) ~= 2),
      %no log file!
      return;
   end
   
   fid = fopen(LogFile,'r');
   if (fid < 0) return, end
   c = fscanf(fid, '%c');
   l = strfind(c,'Can''t bind? tcp_listen[bind]: Address already in use');
   fclose(fid);
   if (~isempty(l)), 
      if (verb) 
         fprintf(2,'Warning: Can''t listen. You probably have another AFNI listening.\nCommunication might fail.\nOnly one AFNI can communicate via plugouts\n');
      end
      return; 
   end
   l = strfind(c,'Plugouts');
   if (~isempty(l)), 
         l = strfind(c,'= listening for connections');
         if (isempty(l)), 
            if (verb) fprintf(2,'Warning: pugouts not enabled. Communication might fail.\n'); end
            return; 
         end
   else 
      %not there yet
      return;      
   end
err = 0;
return;
