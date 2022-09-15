function [err,data] = afni_ni_checkfordata()
   global comm;

   data = [];
   err = 1;
   if (comm.sockcon < 0),
      comm.sockcon = pnet('tcpsocket',sprintf('%d', comm.socknum));
      if (comm.sockcon < 0),
         fprintf(2,'Could not open socket %d\n', comm.socknum);
         return;
      end
      if (comm.dbg ),
         fprintf(1,  'have new socket connection %d\n', comm.sockcon);
      end
   else
      if (comm.dbg > 1) fprintf(1,'have socket %d\n', comm.sockcon); end
   end
   %have socket, try connection
   if (comm.con < 0),
      comm.con = pnet(comm.sockcon,'tcplisten','noblock');  % Don't wait forever
      err = 0;
      if (comm.dbg > 1) fprintf(1,'No connection yet \n'); end
      return;  %nothing yet
   else
      if (comm.dbg > 1) fprintf(1,'have connection handle %d\n', comm.con); end
   end
   status = pnet(comm.con,'status');
   if (status == 0),
      % bad status, perhaps disconnected or first time, reopen socket
      fprintf(2,'Stream closed, attempt to reopen on second pass.\n');
      comm.con = -1;
      err = 0;
      return;
   else
      if (comm.dbg > 1) fprintf(2,'status %d\n', status); end
   end

   %have connection, get data
   pnet(comm.con,'setreadtimeout',0.25);
      % Don't wait forever, I had 3 seconds but that was too slow
      % 0.25 worked just fine on localhost 127.0.0.1

   %check for stuff
   b = pnet(comm.con,'read');
   while (length(b)) data=[data b]; b = pnet(comm.con,'read');  end

   err = 0;
   return
