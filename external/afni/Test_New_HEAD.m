%A script for testing function New_HEAD
%Choose from 3 options

   fprintf(1,'\n Demo for New_HEAD function.\n');
   fprintf(1,'  BEFORE you proceed, make sure you have no AFNI\n');
   fprintf(1,'  session running with the -yesplugouts option.\n');
   fprintf(1,'  1- From a matlab 3 dimensional matrix, create\n');
   fprintf(1,'     afni volumes in short and float formats.\n');
   fprintf(1,'  2- From a matlab 4 dimensional matrix, create\n');
   fprintf(1,'     afni volumes as bucket or timeseries datasets.\n');
   fprintf(1,'  3- From a matlab 3 dimensional matrix, create an\n');
   fprintf(1,'     AFNI volume in Talairach (TLRC) space. The matrix \n');
   fprintf(1,'     data will fill the entire TLRC box.\n');
   iopt = input('\n Enter test number [1 ..3]: ');
   if (iopt ~= 1 & iopt ~= 2 & iopt ~= 3),
      fprintf(1,'Error: must choose either 1, 2 or 3. Not %g\n', iopt);
      return;
   end
   opt = sprintf('test%d', iopt);
   FuncName = 'Test_New_HEAD';
   ppp = sprintf('Test_%d_New_HEAD', iopt);
   if (strcmp(opt,'test1')),
      echo on
      %say you have a matrix of a particular size
      M = flow(100);
      %now to put it in a dataset
      optt.prefix = sprintf('%s', ppp);
      unix(sprintf('rm -f %s*.????  >& /dev/null', optt.prefix));
      optt.dimen  = size(M);
      optt.orient = 'RAI';
      [err,Info1, optt] = New_HEAD(optt); %note that we're also returning optt , some fields get added in New_HEAD for convenience
      if (err),
         fprintf(1,'Error %s:\nFailed in New_HEAD\n', FuncName);
         return;
      end
      [e,m,i] = WriteBrik(M,Info1,optt);
      
      %now write the output as floats
      optt.prefix = sprintf('%s_float', ppp);
      unix(sprintf('rm -f %s*.????  >& /dev/null', optt.prefix));
      optt.datum = 'float';
      [err,Info2, optt] = New_HEAD(optt);
      if (err),
         fprintf(1,'Error %s:\nFailed in New_HEAD\n', FuncName);
         return;
      end

      [e,m,i] = WriteBrik(M,Info2,optt);
      
      echo off
      fprintf(1,'\nAutomatic viewing of results:\n\n');
      cs = []; 
      cs = NewCs('start_afni', '', sprintf('%s.HEAD %s.HEAD ', Info1.RootName, Info2.RootName));
      err = TellAfni(cs); clear cs
      if (err),
         fprintf(2,'Error: Failed to start AFNI in listening mode.\n');
         return;
      end
      i = 1; 
      cs(i) = NewCs('Set_Function', '', Info2.RootName); i = i + 1;
      cs(i) = NewCs('See_Overlay', '', '+'); i = i + 1;
      cs(i) = NewCs('open_window', '', 'axialimage', 'mont=2x2:8 keypress=v geom=500x500+800+50'); i = i+1;
      err = TellAfni(cs); clear cs;
      fprintf(1,'\n\n');
      fprintf(1,'Examine the data interactively.\n');
      fprintf(1,'When done with demo, hit "enter" to quit.\n');
      input (': ','s');
      err = TellAfni(NewCs('Quit'));
      if (err),
         fprintf(2,'Error: Failed telling AFNI.\n');
         return;
      end
      return;
   elseif (strcmp(opt,'test2')),
      echo on
      %Now let us say we have a timeseries of matrices
      n=30; N = 40;
      M = zeros(n,2*n,n,N);
      for (i=1:1:N),
         M(:,:,:,i) = flow(n).*cos(i./N*2.0.*pi)+randn(n,2*n,n)./2.0;
      end
      %Create a multi-brick dataset
      optt.prefix = sprintf('%s_bucket', ppp);
      unix(sprintf('rm -f %s*.????  >& /dev/null', optt.prefix));
      optt.dimen  = size(M);
      optt.orient = 'RAI';
      [err,Info1, optt] = New_HEAD(optt); %note that we're also returning optt , some fields get added in New_HEAD for convenience
      if (err),
         fprintf(1,'Error %s:\nFailed in New_HEAD\n', FuncName);
         return;
      end
      %just put some labels to test label checking 
      Info1.BRICK_LABS = sprintf('%c~', [0:1:Info1.DATASET_RANK(2)-1]+'A');
      Info1.BRICK_LABS = sprintf('''%s', Info1.BRICK_LABS);
      [e,m,i] = WriteBrik(M,Info1,optt);
      %make output a time series of TR 2.0s
      optt.tr = 2.0;
      optt.prefix = sprintf('%s_TS', ppp);
      unix(sprintf('rm -f %s*.????  >& /dev/null', optt.prefix));
      optt.dimen  = size(M);
      optt.orient = 'RAI';
      [err,Info2, optt] = New_HEAD(optt); %note that we're also returning optt , some fields get added in New_HEAD for convenience
      if (err),
         fprintf(1,'Error %s:\nFailed in New_HEAD\n', FuncName);
         return;
      end
      [e,m,i] = WriteBrik(M,Info2,optt);
      echo off
      fprintf(1,'\nAutomatic viewing of results:\n\n');
      cs = []; 
      cs = NewCs('start_afni', '', sprintf('%s.HEAD %s.HEAD ', Info1.RootName, Info2.RootName));
      err = TellAfni(cs); clear cs
      if (err),
         fprintf(2,'Error: Failed to start AFNI in listening mode.\n');
         return;
      end
      i = 1; 
      cs(i) = NewCs('open_window', '', 'axialimage', 'mont=4x4:8 keypress=v geom=500x500+800+50'); i = i+1;
      cs(i) = NewCs('open_window', '', 'axialgraph', 'keypress=v geom=500x500+500+500'); i = i+1;
      err = TellAfni(cs); clear cs;
      fprintf(1,'\n\n');
      fprintf(1,'Examine the data interactively.\n');
      fprintf(1,'When done with demo, hit "enter" to quit.\n');
      input (': ','s');
      err = TellAfni(NewCs('Quit'));
      if (err),
         fprintf(2,'Error: Failed telling AFNI.\n');
         return;
      end
      return;
   elseif (strcmp(opt,'test3')),
      echo on
      %say you have a matrix of a particular size
      M = flow(20);
      %now to put it in a TLRC box 
      optt.view = '+tlrc';       %note that nothing special is done to the coords here
      optt.prefix = sprintf('%s', ppp);
      unix(sprintf('rm -f %s*.????  >& /dev/null', optt.prefix));
      optt.dimen  = size(M);
      [err,Info3, optt] = New_HEAD(optt); %note that we're also returning optt , some fields get added in New_HEAD for convenience
      if (err),
         fprintf(1,'Error %s:\nFailed in New_HEAD\n', FuncName);
         return;
      end
      
      [e,m,i] = WriteBrik(M,Info3,optt);
       
      echo off
      fprintf(1,'\nAutomatic viewing of results:\n\n');
      cs = []; 
      cs = NewCs('start_afni', '', sprintf('%s.HEAD ', Info3.RootName));
      err = TellAfni(cs); clear cs
      if (err),
         fprintf(2,'Error: Failed to start AFNI in listening mode.\n');
         return;
      end
      i = 1; 
      cs(i) = NewCs('open_window', '', 'axialimage', 'geom=500x500+800+50'); i = i+1;
      err = TellAfni(cs); clear cs;
      fprintf(1,'\n\n');
      fprintf(1,'Examine the data interactively.\n');
      fprintf(1,'When done with demo, hit "enter" to quit.\n');
      input (': ','s');
      err = TellAfni(NewCs('Quit'));
      if (err),
         fprintf(2,'Error: Failed telling AFNI.\n');
         return;
      end
      return;
   else
      fprintf(1,'Error %s:\Test = %d is invalid\n', FuncName, iopt);
      return;
   end
return;
