function compile_mex_list(L, baseDir, force)
% function compile_mex_list(L, baseDir)
%
% Compile a list of MEX files as determined by the input argument L.
% The second argument 'baseDir' is the common base directory for the
% files listed in L. The third argument is a flag that determines
% whether to force (re-)compilation even if the MEX file is up-to-date.
%
% See also ft_compile_mex, add_mex_source.

% (C) 2010 S. Klanke

for i=1:length(L)
   [relDir, name] = fileparts(L(i).relName);

   sfname = [baseDir filesep L(i).dir filesep L(i).relName '.c'];
   SF = dir(sfname);
   if numel(SF)<1
      fprintf(1,'Error: source file %s cannot be found.', sfname);
      continue;
   end
   
   if ~force
      mfname = [baseDir filesep L(i).dir filesep name '.' mexext];
      MF = dir(mfname);
      if numel(MF)==1 && datenum(SF.date) <= datenum(MF.date)
         fprintf(1,'Skipping up-to-date MEX file %s/%s\n', L(i).dir, name);
         continue;
      end 
   end
   fprintf(1,'Compiling MEX file %s/%s ...\n', L(i).dir, name);
   cd([baseDir '/' L(i).dir]);
   cmd = sprintf('mex %s.c %s', L(i). relName, L(i).extras);
   eval(cmd);
end
