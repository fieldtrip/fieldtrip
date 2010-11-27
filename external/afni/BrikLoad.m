function [err, V, Info, ErrMessage] = BrikLoad (BrikName, param1, param2)
%
%  OLD USAGE for backward compatibility, see New Usage
%
%   [err, V, Info, ErrMessage] = BrikLoad (BrikName, [form], [MachineFormat])
% or 
%  [V,Info] = BrikLoad (BrikName, [form], [MachineFormat])
% or 
%  [V] = BrikLoad (BrikName, [form], [MachineFormat])
%
%Purpose:
%   loads an AFNI brik into V
%   
%   
%Input Parameters:
%   BrikName, name of the brik
%   form, 'vector' or 'matrix' this is optional, the default is matrix
%        see .Format help under NEW USAGE
%   MachineFormat is a string such as 'native' or 'ieee-le' (LSB_FIRST) or 'ieee-be' (MSB_FIRST)
%        see .MachineFormat help under NEW USAGE
%
%
%  NEW USAGE
%
%   [err, V, Info, ErrMessage] = BrikLoad (BrikName, [Opt])
%
%Purpose:
%   loads an AFNI brik or a NIFTI or a 1D file into V
%   (NIFTI files are loaded via a hidden 3dcopy command)
%   
%Input Parameters:
%   BrikName, name of the brik
%   Opt is an options format with the following optional fields
%   .Format, 'vector' or 'matrix', the default is matrix
%        This determines how the brick data is returned in V.
%        If you choose 'vector' then a N x M x K volume is stored in a 
%        (N * M * K) column vector. If the brick has multiple volumes (4-D)
%        then an N x M x K x J brick is stored in an (N * M * K) x J  matrix.
%        If you use 'vector' option you can change V to matrix format by using
%        M = reshape(...
%              V, Info.DATASET_DIMENSIONS(1), Info.DATASET_DIMENSIONS(2),...
%              Info.DATASET_DIMENSIONS(3), Info.DATASET_RANK(2));
%        
%        Note that indexing in matlab is different than in AFNI. 
%        Matlab starts indexing at  1 while AFNI starts with 0. 
%        So voxel (i,j,k) in AFNI corresponds to voxel (i+1,j+1,k+1) 
%        in matlab. (see example below).
%
%   .MachineFormat is a string such as 'native' or 'ieee-le' (LSB_FIRST) 
%       or 'ieee-be' (MSB_FIRST).
%       Default is whatever is specified in the .HEAD file.
%       If nothing is specified in the .HEAD file, MSB_FIRST or 'ieee-be' 
%       is assumed since most files created with the old versions of AFNI were on %       SGIs. You must specify the parameter Opt.Format, to specify 
%       Opt.MachineFormat and override the defaults.
%       See help fopen for more info 
%   .OutPrecision: If specified, return V with a certain precision. 
%                  The default is to return V in double precision. 
%                  Here are the available options:
%                  '': Default, returns V as double
%                  '*': Returns V in the precision of the Brik itself. 
%          NOTE:        Scaling factors are not applied with this option
%          -----        because of the risk of overflow. You will have to scale
%                       each sub-brick of V by its scaling factor (if any) 
%                       outside of this BrikLoad function. For instructions on 
%                       how to scale a sub-brick, see the section in BrikLoad.m 
%                       under the 'if (Opt.Scale)' statement.
%   .Scale 0/1 if 1, then the scaling factor is applied to the values read 
%        from the brik. Default is 1.
%        Note that Scale cannot be used with OutPrecision.
%     WARNING: If you use .Scale = 0 , you may get bad time series if you have 
%              multiple subbriks ! Each subbrik can have a different scaling 
%              factor
%
%   .Slices: vector of slices, 1 based. Default is all slices. 
%            Read the set of slices specified in .Slices (added for FMRISTAT)
%   .SliceSize_1D: If you are reading 1D files in chunks, specify the 
%                  number of rows you want read in at any one time.
%                  If chunk is 1000 and .Slices = 1 then you would 
%                  get values from rows 1 (the first row) to row 1000.
%                  When .Slices = 2, you would get rows 1001 ... 2000 .
%                  If the number of rows in your dataset is not an 
%                  integer multiple of SliceSize_1D the last read
%                  will return the left over rows.
%   .Frames: vector of frames, 1 based. Default is all frames. 
%            Read the sub-bricks specified in .Frames (added for FMRISTAT)
%   .method: method option for Read_1D if you are using 1D files.
%           see Read_1D -help for more info
%   -------------------------------------------------------------------------
%   (Non-standard options, implemented by Nick Oosterhof (NNO),
%    noosterh@princeton.edu, n.oosterhof@bangor.ac.uk 
%    ** If these options fail, you know who to blame ** )
%
%   .PixX   vector of pixels to read in x direction
%   .PixY   "                         " y "       "
%   Note for .PixX and .PixY: This is intended to reduce the amount of RAM
%   needed when you are only interested in a certain region of the brain.
%   Thus, compared to reading entire volumes, this method may reduce 
%   swapping but will cause more hard drive seeking when reading the data. 
%   In the current implementation, either both fields must be present or 
%   both must be absent. 
%   -------------------------------------------------------------------------
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   V : the brik data in vector or matrix (X,Y,Z,t) form
%   Info : is a structure containing some Brik infor (from .HEAD file)
%     ErrMessage: An error or warning message
%      
%Key Terms:
%   
%More Info :
%   How to read a brick and display a time series:
%   %read the 3D+time brick
%   [err, V, Info, ErrMessage] = BrikLoad ('ARzs_CW_avvr+orig'); 
%   %plot the time course of voxel (29, 33, 3)
%   plot (squeeze(V(30 , 34, 4, :))); 
%   
%   
%
%   see also BrikInfo, Info_1D, WriteBrik
%
%     The core for this function (fread and reshape) is based on 
%     Timothy M. Ellmore's (Laboratory of Brain and Cognition, NIMH)
%     function ReadBRIK 
%     
%     .Slices and .Frames options were added by Dr. Keith Worsley for FMRISTAT
%     http://www.math.mcgill.ca/keith
%     http://www.math.mcgill.ca/keith/fmristat/
%
%     .PixX and .PixY options were added by Nick Oosterhof
%     noosterh@princeton.edu, n.oosterhof@bangor.ac.uk
% 
%     Author : Ziad Saad  Mon Oct 18 14:47:32 CDT 1999 
%     Biomedical Engineering, Marquette University
%     Biophysics Research institute, Medical College of Wisconsin
%
%     Contact: saadz@mail.nih.gov


%Define the function name for easy referencing
FuncName = 'BrikLoad';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
V = [];
Info = [];
ErrMessage = '';

switch nargin,
   case 1,
      DoVect = 0;
      Opt.MachineFormat = '';
      Opt.Format = 'matrix';
      Opt.Scale = 1;
      Opt.Slices = [];
      Opt.Frames = [];
      Opt.FileFormat = '';
      Opt.OutPrecision = '';
   case 2,
      if (~isstruct(param1)),
         Opt.Format = param1;
         Opt.MachineFormat = '';
         Opt.Scale = 1;
         Opt.Slices = [];
         Opt.Frames = [];
         Opt.FileFormat = '';
         Opt.OutPrecision = '';
      else
         Opt = param1;
         if (~isfield(Opt,'FileFormat')),   Opt.FileFormat = ''; end
         if (~isfield(Opt,'OutPrecision')),   Opt.OutPrecision = ''; end
      end
   case 3,
         Opt.Format = param1;
         Opt.MachineFormat = param2;
         Opt.Scale = 1;
         Opt.Slices = [];
         Opt.Frames = []; 
         Opt.FileFormat = '';
         Opt.OutPrecision = '';         
   otherwise,
   err= ErrEval(FuncName,'Err_Bad option');
   return;
end

% Are we dealing with a .1D file ?
is1D = 0;
if (isempty(Opt.FileFormat)), % try to guess 
   [St, xtr] = Remove1DExtension(BrikName);
   if (~isempty(xtr)),
      is1D = 1;
   end
elseif (strcmp(Opt.FileFormat, '1D') | strcmp(Opt.FileFormat, '1d')),
   is1D = 1;
end
isNIFTI = 0;
if (isempty(Opt.FileFormat)), % try to guess 
   [St, xtr] = RemoveNiftiExtension(BrikName);
   if (~isempty(xtr)),
      isNIFTI = 1;
   end
elseif (   strcmp(Opt.FileFormat, 'NIFTI') ...
         | strcmp(Opt.FileFormat, 'nifti') ...
         | strcmp(Opt.FileFormat, 'Nifti') ...
         | strcmp(Opt.FileFormat, 'Nii') ...
         | strcmp(Opt.FileFormat, 'nii') ...
         | strcmp(Opt.FileFormat, 'NII')),
   isNIFTI = 1;
end

if (is1D), % 1D land
   V = []; Info = []; ErrMessage = '';
   Opt.verb = 1;
   if (~isfield(Opt, 'method') | isempty(Opt.method)), Opt.method = 0; end
   if (isfield(Opt,'Slices') & ~isempty(Opt.Slices)),
      if (length(Opt.Slices) ~= 1),
         ErrMessage = sprintf ('%s: Opt.Slices can only be used to specify one slice at a time for 1D format', FuncName, BrikName);
         err = ErrEval(FuncName,'Err_1D Bad .Slices option');
         return;
      end
      if (~isfield(Opt, 'SliceSize_1D') | isempty(Opt.SliceSize_1D)),
         ErrMessage = sprintf ('%s: SliceSize_1D must be specified with Slices option for 1D files', FuncName);
         err = ErrEval(FuncName,'Err_1D Bad .SliceSize_1D option');
         return;
      end
      Opt.chunk_size = Opt.SliceSize_1D;
      Opt.chunk_index = Opt.Slices - 1;
      if (isfield(Opt,'Frames') & ~isempty(Opt.Frames)), Opt.col_index = Opt.Frames - 1; end
   end
   
   [err, V, Info] = Read_1D(BrikName, Opt);
   if (err), 
      ErrMessage = sprintf ('%s: Failed to read %s file', FuncName, BrikName);
      err = ErrEval(FuncName,'Err_1D file could not be read');
      return;
   end
   if (nargout == 4 | nargout == 3),   
      err = 0;
   else 
      err = V;
      V = Info;
   end

   return;
end

if (0), %fails with NIFTI names including paths, 
        %delete this section in near future 
   if (isNIFTI),
      [St, xtr] = RemoveNiftiExtension(BrikName); 
      NiftiPref = St;          
      %create a BRIK version
      otmp = sprintf('./____tmp_%s', St);
      stmp = sprintf('3dcopy %s %s', BrikName, otmp);
      [us, uw] = unix(stmp);
      if (us),
         ErrMessage = sprintf ('%s: Failed to create afni brick:\n%s',...
                               FuncName, uw);
         err = ErrEval(FuncName,'Err_Could not create afni brick');
         return;
      end
      obrik = dir(sprintf('%s+*',otmp));
      if (length(obrik) ~= 2),
         ErrMessage = sprintf ('%s: Failed to find afni brick', FuncName);
         err = ErrEval(FuncName,'Err_Could not find afni brick');
         return;
      end
      BrikName = obrik(1).name;
   end
else, %Craig Stark's fix for cases with directory names 
   if (isNIFTI)
      [St, xtr] = RemoveNiftiExtension(BrikName); 
      NiftiPref = St;          
      %create a BRIK version
      subdir = '';
      if ~isempty(strfind(St,'/')) % We have a directory name in here
         dir_pos = strfind(St,'/');
         dir_pos = dir_pos(length(dir_pos)); % Get the last /
         subdir = St(1:dir_pos);
         otmp = sprintf('%s____tmp_%s',subdir,St((dir_pos+1):length(St)));
      else
         otmp = sprintf('./____tmp_%s', St);
      end
      stmp = sprintf('3dcopy %s %s', BrikName, otmp);
      [us, uw] = unix(stmp);
      if (us),
         ErrMessage = sprintf ('%s: Failed to create afni brick:\n%s',...
                               FuncName, uw);
         err = ErrEval(FuncName,'Err_Could not create afni brick');
         return;
      end
      obrik = dir(sprintf('%s+*',otmp));
      if (length(obrik) ~= 2),
         ErrMessage = sprintf ('%s: Failed to find afni brick', FuncName);
         err = ErrEval(FuncName,'Err_Could not find afni brick');
         return;
      end
      obrik(1).name = sprintf('%s%s',subdir,obrik(1).name);
      obrik(2).name = sprintf('%s%s',subdir,obrik(2).name);
      BrikName = obrik(1).name;
   end
end

%assume you have a brik format and proceed as usual 
%make sure Opt fields are set OK. That's done for the new usage format
%   if (~isfield(Opt,'') | isempty(Opt.)),   Opt. = ; end
   if (~isfield(Opt,'Format') | isempty(Opt.Format)),   Opt.Format = 'matrix'; end
   if (~isfield(Opt,'MachineFormat') | isempty(Opt.MachineFormat)),   Opt.MachineFormat = ''; end
   if (~isfield(Opt,'Scale') | isempty(Opt.Scale)),   Opt.Scale = 1; end %you can change the default for the new usage here
   if (~isfield(Opt,'Slices') | isempty(Opt.Slices)),   Opt.Slices = []; end 
   if (~isfield(Opt,'Frames') | isempty(Opt.Frames)),   Opt.Frames = []; end 
   if (~isfield(Opt,'FileFormat') | isempty(Opt.FileFormat)),   Opt.FileFormat = ''; end
   if (~isfield(Opt,'OutPrecision') | isempty(Opt.OutPrecision)),   Opt.OutPrecision = ''; end

%Check for conflicts between OutPrecision and Scale
if (Opt.Scale & ~isempty(Opt.OutPrecision)),
   ErrMessage = sprintf ('%s: Opt.Scale = 1 cannot be used with Opt.OutPrecision\n', FuncName);
   err = ErrEval(FuncName,'Err_Opt.Scale = 1 cannot be used with Opt.OutPrecision');
   return;
end

%check on OutPrecision
if (~isempty(Opt.OutPrecision)),
   if (~strcmp(Opt.OutPrecision,'*')),
      ErrMessage = sprintf ('%s: Allowed OutPrecision value is ''*''\n', FuncName);
      err = ErrEval(FuncName,'Err_Bad value for Opt.OutPrecision');
   end
end

%set the format   
   if (eq_str(Opt.Format,'vector')),
      DoVect = 1;
   elseif(eq_str(Opt.Format,'matrix')),
      DoVect = 0;
   end

%Fix the name of the brik
   %make sure there are no .BRIK or .HEAD   
   vtmp = findstr(BrikName,'.BRIK');
   if (~isempty(vtmp)), %remove .BRIK
      BrikName = BrikName(1:vtmp(1)-1);
   end
   vtmp = findstr(BrikName,'.HEAD');
   if (~isempty(vtmp)), %remove .HEAD
      BrikName = BrikName(1:vtmp(1)-1);
   end
   
   %remove last dot if it is in name
   if (BrikName(length(BrikName)) == '.' & length(BrikName)>1),
      BrikName = BrikName(1:length(BrikName)-1);
   end

   %Now make sure .HEAD and .BRIK are present
   sHead = sprintf('%s.HEAD', BrikName);
   sBRIK = sprintf('%s.BRIK', BrikName);
   if (~filexist(sHead)), 
      ErrMessage = sprintf ('%s: %s not found', FuncName, sHead);
      err = ErrEval(FuncName,sprintf ('Err_HEAD file %s not found', sHead));
      return;
   end
   sBRIK_gz = sprintf('%s.BRIK.gz', BrikName);   
   if (filexist(sBRIK_gz)),
      fprintf(2,'Unzipping %s...\n', sBRIK_gz);
      unix(sprintf('gzip -d %s', sBRIK_gz)); %matlab gunzip version requires Java
   end
   if (~filexist(sBRIK)), 
      ErrMessage = sprintf ('%s: %s not found', FuncName, sBRIK);
      err = ErrEval(FuncName,'Err_BRIK file not found');
      return;
   end

   
%get Brik info and setup for slice and frame selectors (for FMRISTAT)
   [err, Info] = BrikInfo(BrikName);

   allslices=1:Info.DATASET_DIMENSIONS(3);
   if  (~isfield(Opt,'Slices') | isempty(Opt.Slices)),   
      Opt.Slices = allslices; 
   end
   isallslices=all(ismember(allslices,Opt.Slices));
   allframes=1:Info.DATASET_RANK(2);
   if  (~isfield(Opt,'Frames') | isempty(Opt.Frames)),   
      Opt.Frames = allframes; 
   end
   isallframes=all(ismember(allframes,Opt.Frames));

   if (~strcmp(Info.ByteOrder,'unspecified')),
      %found byte order specs, make sure it does not conflict with the user option
      if (~isempty(Opt.MachineFormat)),
         %user specified a format
         if (~strcmp(Info.ByteOrder,Opt.MachineFormat)),
            %clash, warn
            ErrEval(FuncName,'Wrn_Machine format specified conflicts with .HEAD info. Proceeding with fread ...');
         end
      else
         Opt.MachineFormat = Info.ByteOrder;
      end
   else
      %did not find ByteOrder in .HEAD, use user option or pick default
      if (isempty(Opt.MachineFormat)),
         Opt.MachineFormat = 'ieee-be';
         ErrEval(FuncName,'Wrn_No Machine Format was specified by user or found in .HEAD file. Using ieee-be (MSB_FIRST)');
      end
   end

%figure out storage type
   itype = unique(Info.BRICK_TYPES);
   if (length(itype) > 1),
      err =  1; ErrMessage = sprintf('Error %s: Not set up to read sub-bricks of multiple sub-types', FuncName); errordlg(ErrMessage);
      return;
   end
   switch itype,
      case 0
         typestr = 'ubit8';
      case 1
         typestr = 'short';
      case 2
         typestr = 'int';
      case 3
         typestr = 'float';
      otherwise
         err = ErrEval(FuncName,'Err_Cannot read this data type');
         return;
   end
   
%read .BRIK file
   BrikName = sprintf('%s.BRIK', BrikName);
   
   fidBRIK = fopen(BrikName, 'rb', Opt.MachineFormat);
   if (fidBRIK < 0),
      err = ErrEval(FuncName,'Err_Could not open .BRIK file');
      return;
   end
   
   %NNO added
   haspixX=isfield(Opt, 'PixX') && ~isempty(Opt.PixX);
   haspixY=isfield(Opt, 'PixY') && ~isempty(Opt.PixY);
   
   if xor(haspixX, haspixY) %either both or neither
       err =  1; ErrMessage = sprintf('Error %s: options PixX and PixY must be either both present or both absent ', FuncName); errordlg(ErrMessage);
      return;
   end
   
   allpixels=~haspixX && ~haspixY; %(the && is not necessary)
   if allpixels
       numpixX=Info.DATASET_DIMENSIONS(1);
       numpixY=Info.DATASET_DIMENSIONS(2);
       Opt.PixX = [1:1:Info.DATASET_DIMENSIONS(1)];
       Opt.PixY = [1:1:Info.DATASET_DIMENSIONS(2)];
   else
       numpixX=length(Opt.PixX);
       numpixY=length(Opt.PixY);
   end
   numpix=numpixX*numpixY;
   
   numslices=length(Opt.Slices);
   numframes=length(Opt.Frames);

   if isallslices && isallframes && allpixels
      V = fread(fidBRIK, (Info.DATASET_DIMENSIONS(1) .* Info.DATASET_DIMENSIONS(2) .* Info.DATASET_DIMENSIONS(3) .* Info.DATASET_RANK(2)) , [Opt.OutPrecision,typestr]);
   else
      V=zeros(1,numpix*numslices*numframes);
      for k=1:numframes
         frame=Opt.Frames(k);
         if isallslices && allpixels
            fseek(fidBRIK, numpix*Info.DATASET_DIMENSIONS(3)*(frame-1)*Info.TypeBytes, 'bof');
            istrt=1+numpix*numslices*(k-1);
            istp=istrt-1+numpix*numslices;
            V(istrt:istp)=fread(fidBRIK,numpix*numslices,[Opt.OutPrecision,typestr]);
         else
            consecutivex=isequal(min(Opt.PixX):max(Opt.PixX), Opt.PixX); %true iff Opt.PixX = [n, n+1, ..., n+k] 
            
            for j=1:numslices
                slice=Opt.Slices(j);
                if allpixels %NNO added
                    seekidx=numpix*(slice-1+Info.DATASET_DIMENSIONS(3)*(frame-1));
                    fseek(fidBRIK, seekidx*Info.TypeBytes, 'bof');
                    istrt=1+numpix*(j-1+numslices*(k-1));
                    istp=istrt-1+numpix;
                    V(istrt:istp)=fread(fidBRIK,numpix,[Opt.OutPrecision,typestr]);
                else %NNO - read pixels within the slice
                    dims=Info.DATASET_DIMENSIONS;

                    for q=1:numpixY
                        y=Opt.PixY(q);
                        
                        if consecutivex %read all voxels in this row in one go
                            x=Opt.PixX(1);
                            
                            seekidx=(y-1)*dims(1)+(x-1)+dims(1)*dims(2)*(slice-1+dims(3)*(frame-1));
                            fseek(fidBRIK, seekidx*Info.TypeBytes, 'bof');
                            
                            data=fread(fidBRIK,numpixX,[Opt.OutPrecision,typestr]);
                            istrt=numpix*(j-1+numslices*(k-1))+(q-1)*numpixX+1;
                            istp=istrt+numpixX-1;
                            V(istrt:istp)=data;
                        else %the very slowest method possible - read byte by byte
                            for p=1:numpixX
                                x=Opt.PixX(p);

                                seekidx=(y-1)*dims(1)+(x-1)+dims(1)*dims(2)*(slice-1+dims(3)*(frame-1));
                                fseek(fidBRIK, seekidx*Info.TypeBytes, 'bof');

                                data=fread(fidBRIK,1,[Opt.OutPrecision,typestr]);
                                istrt=numpix*(j-1+numslices*(k-1))+(q-1)*numpixX+p;
                                V(istrt)=data;
                            end
                        end
                    end
                end
            end
         end
      end
   end
   fclose (fidBRIK);

%scale if required
   if (Opt.Scale),
      facs=Info.BRICK_FLOAT_FACS(Opt.Frames);
      iscl = find (facs); %find non zero scales
      NperBrik = numpix * numslices; %NNO Used DATASET_DIMENSIONS(1 and 2)
      if (~isempty(iscl)),
         for j=1:1:length(iscl),
            istrt = 1+ (iscl(j)-1).*NperBrik;
            istp = istrt + NperBrik - 1;
            V(istrt:istp) = V(istrt:istp) .* facs(iscl(j));
         end
      end
   else, %give warning
      facs=Info.BRICK_FLOAT_FACS(Opt.Frames);
      iscl = find (facs); %find non zero scales
      if (~isempty(iscl)),
         fprintf(2,'WARNING: Some sub-bricks read have non-zero scaling factors that were not applied\n');
         fprintf(2,'         because Opt.Scale is turned off.\n');
         fprintf(2,'         Those scale factors are in Info.BRICK_FLOAT_FACS and likely should be\n');
         fprintf(2,'         applied before using values in V\n');
      end
   end

%NNO
if DoVect,
    V = reshape(V, numpixX * numpixY * numslices, numframes);
else
    V = reshape(V, numpixX, numpixY, numslices, numframes);
end
    
    
if (isNIFTI),
   if (filexist(obrik(1).name)) delete(obrik(1).name);
   else
      tmpname = RemoveExtension(obrik(1).name, '.gz');
      delete(tmpname);
   end
   delete(obrik(2).name);
   Info.RootName = NiftiPref;
end

if (nargout == 4 | nargout == 3),   
   err = 0;
else 
   err = V;
   V = Info;
end

return;

%test code
clear all
Opt.Slices=[2 8];
Opt.Frames=[1 2];

BrikName='./r1_time+orig.BRIK';
[err, V, Info, ErrMessage] = BrikLoad(BrikName);
size(V)
figure(1);
k=0;
for i=1:2
   for j=1:2
      k=k+1;
      subplot(2,2,k);
imagesc(V(:,:,Opt.Slices(i),Opt.Frames(j))); 
if (exist('spectral')) colormap(spectral); end
colorbar;
end
end

clear all
Opt.Slices=[2 8];
Opt.Frames=[1 2];


BrikName='./r1_time+orig.BRIK';
[err, V, Info, ErrMessage] = BrikLoad(BrikName, Opt);
figure(2);
k=0;
for i=1:2
   for j=1:2
      k=k+1;
      subplot(2,2,k);
imagesc(V(:,:,i,j)); colormap(spectral); colorbar;
end
end


