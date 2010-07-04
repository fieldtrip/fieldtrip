function [nel, comm] = ToyTalk()
%           ToyTalk()
%
% A sample function to illustrate how matlab can listen to communication from
% SUMA. 
%
%Dependency:
%-----------
%The TCP communication is based on pnet.m by:  
%       Peter Rydesäter, 
%       Mitthögskolan(Mid Sweden University) campus Östersund, SWEDEN
%pnet is not in AFNI's matlab library yet. You should have a copy of it by now
%since you pointed me to it. If the mex file is not available for your machine,
%you will need to compile a mex file for pnet. On the mac (or linux), just do:
%   mex -O pnet.c
%
%Of course you will need the AFNI_matlab library that is included in 
%afni's source archive. 
%
%
%Usage:
%------
%  Prepare some data, by creating a coarse (just for simplicity) version
%  of the surfaces used in the AFNI workshop. 
%
%  MapIcosahedron -spec DemoSubj_lh.spec -ld 12 -morph sphere.reg -prefix std12.
%
%  Launch SUMA from the command line using something like:
%
%     setenv SUMA_NI_TEXT_TALK_MODE YES
%
%     suma -i_fs std12.lh.smoothwm.asc -sv DemoSubj_SurfVol+orig. -ah 127.0.0.1
%
%       Note that for now, you need both the -sv AND the -ah options.
%        In future releases this will not be the case. Any volume for -sv would
%        work, if none are relevant for the surface you're using.
%
%  In matlab, run: ToyTalk()
%
%  In SUMA, press 't'
%
%  In matlab, watch it deal with incoming information.
%
%  In SUMA, click around
%
%  In matlab, watch it process the cross hair information.
%
%  To stop the function, and get the matlab prompt again, 
%  use ctrl+c 'in matlab'. You can start it again, but will have
%  to press 't' in the still running SUMA twice. 
%  
%  I am hoping the code is sufficient to show you how this is done
%  I tried to hide all unecessary detail for you. It bugs me that 
%  there is no way to start a workprocess on command line, so that
%  one's code does not have to handle repeated polling for new data
%  but oh well. If you new of new matlab developments on that front, 
%  please let me know.
% 
%Limitations:
%-------------
%1- Matlab uses the same port as AFNI for listening to SUMA. So you should not be
%running AFNI with -niml (listening mode) if you want matlab to be listening too.
%This limitation will soon be removed.
%2- Matlab does not talk 'directly' to SUMA the way afni does. While this is a 
%possiblity, I need concrete examples of where that is useful. I would be more
%inclined to create a 'TellSuma.m' which parallels 'TellAfni.m' to communicate
%with SUMA from matlab. Let me know what you think
%
  
   global comm ;
  
   comm = afni_talk_defs(); %reset communication structure. 
   comm.dbg = 1;
   pnet ('closeall');
   
   pause on;
   
   err = 0;
   dataread=[];
   ud = struct;    % user data structure
   while (~err ),
      %do something
      figure(1); plot (sin(randn(1,100))); drawnow;
      
      %%%%%%%%%%%%%%% THIS BLOCK SHOULD BE A WORKPROCESS 
      %%%%%%%%%%%%%%% that calls a callback with the data
      %%%%%%%%%%%%%%% that was read. But I don't know how
      %%%%%%%%%%%%%%% to do such a thing yet... 
                      
      %check for meat
      [err,dataread] = afni_ni_checkfordata();
      if (length(dataread)),
         expr = '<(\w+).*?>.*?</\1>';
         [ud.tt ud.mm] = regexp(dataread, expr,'tokens', 'match');
         fprintf(1,[ 'Total data: %d characters\n',...
                     '%d elements\n'],...
                      length(dataread), length(ud.tt));
         for (iel=1:1:length(ud.tt)), 
            fprintf(1,'  %d: %s in mm(%d)\n', iel, char(ud.tt{iel}), iel); 
         end
         clear (dataread);
         
         %change into nel elements 
               %NOTE: You cannot assemble output of afni_nel_parse
               %into a vector of structures because each nel
               %is a struct with different fields
         N_nel = length(ud.mm);
         data(N_nel) = struct;
         for (iel=1:1:N_nel),
            data(iel).nel = afni_nel_parse(ud.mm{iel});
         end
         clear ud
         
         %do the callback
         ProcessSUMAdata(data);
         
         %cleanup
         clear data      
      else
         %no data read
      end
            %%%%%%%%%%%%%%% END of WORKPROCESS block
   end
   return
   
function ProcessSUMAdata(ds)
   global SO;  %surface object
   
   dopt.verbose = 0;
   dopt.OpenGL = 1;
   dopt.GraphType = 'Surf';
   
   if (isempty(SO)),
      SO.NodeNormList = [];
      SO.FaceSetList = [];
      SO.NodeList = [];
   end
   
   N_nel = length(ds);
   for (iel = 1:1:N_nel),     %Process each element. Cleanup is handled in
                              %parent function
      fprintf(1,'Got nel %s\n',ds(iel).nel.name); 
      %do something with nels
      if (strcmp(ds(iel).nel.name,'SUMA_crosshair_xyz')),
         %Get the node of interest
         inode = str2num(ds(iel).nel.surface_nodeid);
         figure(2); hold on
         if (inode >= 0 & isfield(SO,'NodeList') & ~isempty(SO.NodeList)),
            xyz = SO.NodeList(inode+1,:);
            plot3(xyz(1), xyz(2), xyz(3), 'b*', 'MarkerSize', 20);
         else
            xyz = ds(iel).nel.data;
            plot3(xyz(1), xyz(2), xyz(3), 'go', 'MarkerSize', 20);
         end
      elseif (strcmp(ds(iel).nel.name,'SUMA_node_normals')), %node normals
         SO.NodeNormList = ds(iel).nel.data(:,1:3);
      elseif (strcmp(ds(iel).nel.name,'SUMA_ixyz')), %node coords
         SO.NodeList = ds(iel).nel.data(:,2:4);
         if (isfield(SO,'FaceSetList') & ~isempty(SO.FaceSetList)),
            figure(2);hold on
            DispIVSurf(SO.NodeList, SO.FaceSetList+1, [], 2, 2, dopt);
         end
      elseif (strcmp(ds(iel).nel.name,'SUMA_ijk')),   %triangulation
         SO.FaceSetList = ds(iel).nel.data(:,1:3);
         if (isfield(SO,'NodeList') & ~isempty(SO.NodeList)),
            figure(2);hold on
            DispIVSurf(SO.NodeList, SO.FaceSetList+1, [], 2, 2, dopt);
         end
      else
         fprintf(1,'Received %s, don''t know what to do with it\n', ...
                  ds(iel).nel.name);
      end
   end
   return;

