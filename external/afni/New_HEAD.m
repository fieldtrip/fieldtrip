function [err,Info, opt] = New_HEAD (opt)
%
%   [err,Info, Opt] = New_HEAD (Opt)
%
%Purpose:
%   A function for the likes of Greg Detre.
%   It makes the creation of an AFNI dataset
%   from a 3 or 4 dimensional matrix in matlab a breeze.
%   
%   
%   
%Input Parameters:
%   opt is an options structure with the following fields:
%     .prefix: A string containg the AFNI dataset's prefix
%              .prefix and one of .dimen or .master must 
%              be used.
%     .dimen: A string or array containing the dimensions of
%             the dataset (at least 3 dimensions need be present).
%     The following four options are from 3dUndump. See 3dUndump -help   
%     .view:    '+orig' (default) or '+acpc' or '+tlrc' 
%               If you choose +tlrc, the resultant volume
%               will fit the Talairach box.
%     .master: Name of an existing AFNI dataset that would
%              provide the needed matrix size, orientation etc.
%              Do not mix .master with .dimen
%     .datum:  'byte', 'short' (default), or 'float'
%     .orient: Orientation of data in matrix. Default in 'RAI'
% 
%     .tr: a float specifying the TR in seconds.
%          The presence of such a field necessitates a 4 dimensional
%     .scale: 0/1 This is only meaningful with 'short' data where it
%             defaults to 1. It defaults to 0 for other types. It is
%             used when storing float data as shorts to minimize disk
%             use while preserving numeric percision.      
%     .Overwrite: y/[n] allow header to be created even if one with
%                 similar name is found on disk
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   Info: The header structure, given the options specified
%   Opt: A modified version of the input Opt. It will be passed
%        along with the data array to WriteBrik for writing the
%        dataset to disk. (See Examples below)
%      
%Examples: 
%     Run New_HEAD('test') or see Test_New_HEAD.m for examples. 
%   
%More Info :
%     AFNI programs 3drefit and 3dAttribute are your friends
%     BrikLoad, BrikInfo, WriteBrik
%     README.attributes
%
%     The function makes heavy use of 3dUndump
%
%This function requires AFNI binaries postdating Feb. 12 2007
%
%
%     Author : Ziad Saad   saadz@mail.nih.gov
%     Date : Fri Feb 9 16:29:51 EST 2007
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'New_HEAD';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Info = [];

%check on return parameters
if (nargout == 2),
   fprintf(1, 'Notice %s:\nYou are not specifying a return argument for Opt.\n', FuncName);
   fprintf(1, ' Be aware that some fields in Opt might have been modified inside\n');
   fprintf(1, ' %s but these changes are not reflected in the Opt structure you\n', FuncName);
   fprintf(1, ' will be passing to WriteBrik.\n');
   fprintf(1, ' If this makes no sense to you, just use:\n');
   fprintf(1, '    [err, Info, Opt] = New_HEAD(Opt);\n'); 
end

%work options
if (ischar(opt)), %perhaps test mode
   Test_New_HEAD;
   return;
end

%standardize input
if (isfield(opt,'dimen')),
   if (ischar(opt.dimen)),
      v = str2num(opt.dimen);
      rmfield(opt, 'dimen');
      opt.dimen = v;
   end
end


%create command for -master option or non-tlrc views
if (~isfield(opt,'view') | ~strcmp(opt.view,'+tlrc') | isfield(opt,'master')),
   tmp_suf = '___NeW_hEaD_';
   sopt = '3dUndump -head_only';
   if (isfield(opt,'dimen')),
      sopt = sprintf('%s -dimen %s', sopt, num2str(opt.dimen(1:3)));
   end
   if (isfield(opt,'orient')),
      sopt = sprintf('%s -orient %s', sopt, opt.orient);
   end
   if (isfield(opt,'master')),
      [Status, mPrefix, mView] = PrefixStatus (opt.master);
      sopt = sprintf('%s -master %s', sopt, opt.master);
   else
      mView = '+orig';
   end
   if (~isfield(opt,'datum')), opt.datum = 'short'; end

   if (isfield(opt,'datum')),
      if (~strcmp(opt.datum,'short') & ~strcmp(opt.datum, 'byte') & ~strcmp(opt.datum,'float')),
         fprintf(1,'Error %s:\ndatum option must be one of ''byte'', ''short'', ''float''\nI have %s\n', FuncName, opt.datum);
         return;
      end
      sopt = sprintf('%s -datum %s', sopt, opt.datum);
   end

   if (~isfield(opt,'view')),
      opt.view = '+orig';
   end

   if (isfield(opt,'prefix')),
      [Status, Prefix, View] = PrefixStatus (opt.prefix);
      if (Status < 1 & strcmp(View,opt.view)),
         if (isfield(opt,'overwrite') & ~strcmp(opt.overwrite,'y')),
            fprintf(1,'Error %s:\nLooks like %s exists already.\n',...
                     FuncName, opt.prefix);
            return;
         end
      end
      opt.prefix = Prefix;
      ohead = sprintf('%s%s', tmp_suf, opt.prefix);
      sopt = sprintf('%s -prefix %s%s', sopt, ohead);
   else
      fprintf(1,'Error %s:\nNeed a .prefix option.\n', FuncName);
      return;
   end

   [e,w] = unix(sopt);
   if (e),
      fprintf(1,'Error %s:\nHeader creating command %s failed.\nSee this function''s help and 3dUndump -help\n3dUndump''s output was:\n%s\n', ...
               FuncName, sopt, w); 
      New_HEAD_CLEAN(tmp_suf); 
      return;  
   end

   [err,Info]  = BrikInfo(sprintf('%s%s', ohead,mView));
   New_HEAD_CLEAN(tmp_suf); 
else %have +tlrc and no master  
   mView = '+tlrc';
   [e,d] = unix('which afni');
   if (e),
      fprintf(1,'Error %s:\nFailed to find afni!\n', FuncName);
      return;  
   end
   [e,pt] = GetPath(d);
   templ = sprintf('%s/TT_N27+tlrc.HEAD', pt);
   [e,Info] = BrikInfo(templ);
   if (e),
      fprintf(1,'Error %s:\nFailed to get header of %s!\n', FuncName, templ);
      return;  
   end
   %recalculate the origin (res. of template is 1x1x1mm)
   oOrig = Info.ORIGIN;
   oDelta = Info.DELTA;
   oDimen = Info.DATASET_DIMENSIONS;
   %Calculate dimension ratio and adjust so that volume fits in box.
   rat = oDimen(1:3)./opt.dimen(1:3);
   Info.DELTA = Info.DELTA .* rat;
   Info.DATASET_DIMENSIONS(1:3) = opt.dimen(1:3);
   %Now shift the origin so that the final volume still fits in the same box (edge to edge) of TLRC box
   dOrig = 1./2*(Info.DELTA-oDelta);
   Info.ORIGIN = Info.ORIGIN+dOrig;
   %Have to deal with datum
   if (isfield(opt,'datum')),
      if (strcmp(opt.datum,'byte')),
         Info.BRICK_TYPES(1) = 0;
      elseif (strcmp(opt.datum,'short')),
         Info.BRICK_TYPES(1) = 1;
      elseif (strcmp(opt.datum,'float')),
         Info.BRICK_TYPES(1) = 3;
      else
         fprintf(1,'Error %s: Bad data type %s\n', FuncName, opt.datum);
         return;
      end
   else
      opt.datum = 'short'
      Info.BRICK_TYPES(1) = 1;
   end
end

%the scaling option
if (isfield(opt,'scale')),
   if (opt.scale & ~strcmp(opt.datum,'short')),
      fprintf(1,'Warning %s:\n .scale option is only for ''short'' type data.\n Resetting it to 0.\n', FuncName);
      opt.scale = 0;
   end
else 
   if (strcmp(opt.datum,'short')), 
      opt.scale = 1;
   else
      opt.scale = 0;
   end
end



%take care of prefix business
Info.RootName = sprintf('%s%s', opt.prefix, opt.view);

%take care of view
if (strcmp(mView, opt.view) == 0), %different
   if (strcmp(opt.view,'+orig')) Info.SCENE_DATA(1) = 0;
   elseif (strcmp(opt.view,'+acpc')) Info.SCENE_DATA(1) = 1;
   elseif (strcmp(opt.view,'+tlrc')) Info.SCENE_DATA(1) = 2;
   else
      fprintf(1,'Error %s:\nBad view %s\n', FuncName, opt.view);
      return;
   end 
end

%take care of 4th dimen
if (isfield(opt,'dimen')),
   if (length(opt.dimen) > 3),
      Info.DATASET_RANK(2) = opt.dimen(4);
      Info.BRICK_TYPES = Info.BRICK_TYPES(1).*ones(1,opt.dimen(4));
      Info.BRICK_FLOAT_FACS  = Info.BRICK_FLOAT_FACS(1) .*ones(1, opt.dimen(4));
   end
end

%take care of TR
if (isfield(opt,'tr')),
   Info.TAXIS_NUMS(1) = Info.DATASET_RANK(2);
   Info.TAXIS_NUMS(2) = 0; %no time offset for slices at the moment
   Info.TAXIS_NUMS(3) = 77002; %units in seconds for tr (below)

   Info.TAXIS_FLOATS(1) = 0;  %time origin 0
   Info.TAXIS_FLOATS(2) = opt.tr;   %TR in units of Info.TAXIS_NUMS(3)
   Info.TAXIS_FLOATS(3) = 0; %duration of acquisition
   Info.TAXIS_FLOATS(4) = 0; %no time offset please
   Info.TAXIS_FLOATS(5) = 0; %no time offset please
   
   Info.TAXIS_OFFSETS = zeros(1,Info.TAXIS_NUMS(1));  %no bloody time offset
end

err = 0;
return;

function New_HEAD_CLEAN(sss)
   unix(sprintf('rm -f %s*.HEAD >& /dev/null', sss));
return;

