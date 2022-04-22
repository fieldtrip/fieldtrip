% comm = afni_talk_defs() for getting communications defaults
% If you have a port offset PO to pass along, use:
%     comm = afni_talk_defs(comm, PO);
function [comm] = afni_talk_defs(comm, np)

   comm.con = -1;
   comm.sockcon = -1;
   if (nargin > 1),
      com = sprintf('suma  -npq %d -port_number_quiet MATLAB_SUMA_NIML', np);
   else
      com = sprintf('suma  -port_number_quiet MATLAB_SUMA_NIML');
   end
   [l,a] = system(com);
   if (l),
      fprintf(2,...
             'Error fetching port number. Assuming default of 53211\n%s\n', a);
      comm.socknum = 53211;
   else
      comm.socknum = str2double(a);
   end
   comm.dbg = 0;
   comm.mm = {};
   comm.tt = {};

   return;
