function [MapMPM, MapML, aur] = CA_EZ_Prep_genx()
% A function to process a new Eickhoff, Amunts and Zilles CA toolbox
% create NIML files that describe "atlas point lists" for each
% of the atlases included
% DRG/ZSS SSCC MAY 2011

FuncName = 'CA_EZ_Prep_genx';
MapMPM = [];
MapML = [];

% version 1.8
CA_EZ_Version = '18'
toolbox_dir = '/Volumes/Data/atlas/eickhoff/Anatomy';
% toolbox_dir = '/Volumes/afni/home4/users/ziad/Programs/matlab/spm2/toolbox/Anatomy';
nimlout_dir = toolbox_dir;

if (exist(toolbox_dir) ~= 7),
   fprintf(2,'Anatomy toolbox directory %s not found\nPick a new one:', toolbox_dir);
   toolbox = uigetdir(sprintf('.%c', filesep), 'Standard toobox dir not found. Pick a new one:');
   if (exist(toolbox_dir) ~= 7),
      fprintf(2,'Toolbox directory %s not found.', toolbox_dir);
      return;
   end
else
   %get around the symbolic linc so that reference is to actual directory ...
   curdir = pwd;
   cd (toolbox_dir);
   toolbox_dir = pwd;
%   cd (curdir);
   [err,pt] = GetPath(toolbox_dir);
   [err,sout] = unix(sprintf('ls -l %s', pt));
   c = input(sprintf('Found toolbox here: %s\nDirectory Listing:\n%s\nEnter "y" to use it, anything else to quit.\n', toolbox_dir, sout),'s');
   if (isempty(c) || ( c(1) ~= 'y' && c(1) ~= 'Y' ) ),
      return;
   end
end

   fprintf(1,'Using toolbox directory %s...\n', toolbox_dir);

%First get the MPM info
   prf = sprintf('%s%cAllAreas_v%s_MPM.mat', toolbox_dir, filesep, CA_EZ_Version);
   MPM_file = prf;
%    [err, ErrMessage, MPM_file] = zglobb ({prf});
%    if (size(MPM_file,1) ~= 1),
%       fprintf(2,'Could not find unique MPM map list in %s\n', toolbox_dir);
%       for (i=1:1:size(MPM_file,1)), MPM_file(i), end
%       return;
%    end

   %load the MPM map structure
   MapMPM = load(prf);
   MapMPM = MapMPM.MAP;

   %checks
   if (~isstruct(MapMPM)),
      fprintf(2,'MapMPM is not a struct\n');
      return;
   end

   fld = 'name';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'GV';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'ref';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'smoothed';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'VOL';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'MaxMap';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'XYZ';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'XYZmm';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'orient';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'Z';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'LR';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'allXYZ';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'allZ';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end
   fld = 'allLR';
   if (~isfield (MapMPM, fld)), fprintf(2,'%s field is not in MapMPM\n', fld); return; end

%Now the MacroLabels
   prf = sprintf('%s%cMacro.mat', toolbox_dir, filesep);
   [err, ErrMessage, ML_file] = zglobb ({prf});
   if (size(ML_file,1) ~= 1),
      fprintf(2,'Could not find unique ML map list in %s\n', toolbox_dir);
      for (i=1:1:size(ML_file,1)), ML_file(i), end
      return;
   end


   %load the MacroLabels
   MapML = load(ML_file(1).name);
   MapML = MapML.Labels;
   if (~iscellstr(MapML)),
      fprintf(2,'MacroLabels variable not the expected cellstr\n');
      return;
   end
   for (i=1:1:length(MapML)),
      MapML(i) = cellstr(fix_string(char(MapML(i))));
   end
   MapML = char(MapML);

%Output files
   MPMname = 'mpm.niml';
   PMname = 'pm.niml';
   MLname = 'ml.niml';
   LRname = 'lr.niml';
   rname = 'refnames.txt';

   NLbl_ML = size(MapML,1);
   MaxLbl_ML = size(MapML,2)+3;
   NLbl_MPM = length(MapMPM);
   MaxLbl_MPM = 0;
   for (i=1:1:NLbl_MPM),
      if (MaxLbl_MPM < length(MapMPM(i).name)), MaxLbl_MPM = length(MapMPM(i).name); end
   end
   MaxLbl_MPM = MaxLbl_MPM+3;
   if (MaxLbl_MPM > 64),
      fprintf(2,'Error: Labels longer than ATLAS_CMAX defined in AFNI src code.\nIncrease limit here and in thd_ttaltas_query.h\n');
      return;
   end

%Do the references
% a horrible mess of an if block below. One hopes for a better solution someday
% still a mess, but don't need to worry about getting info into C code
% just make text file that will be imported into a atlas NIML table file
% The current directory should now be the toolbox directory, so the correct
% se_note.m should be used.
if (~isempty(which('se_note'))),
   fprintf(1,'\nNow trying to get at references using se_note\n');
   se_note;
   k = 0;
   vers = '';
   if (exist('fg')),
      h = get(fg, 'Children');
      cs = [];
      nref = 0;
      l = 1;
      for (i=1:1:length(h)),
         tmp = get(h(i),'String');
         if (~isempty(tmp)),
            k = k + 1;
            cs(k).s = tmp;
            if (iscellstr(cs(k).s) && length(cs(k).s) > 1),
               if (length(cs(k).s) > 5), %papers
                  nref = nref + 1;
                  ref(nref) = k;
                  cs(k).typ = 1;
               else
                  cs(k).typ = -1; %the authors's info
                  au = cellstr(cs(k).s);
               end
            else
               cs(k).typ = 0; %other strings
               ot(l) = cellstr(cs(k).s);
               if (~isempty(strfind(char(ot(l)), 'Version'))),
                  vers = zdeblank(char(ot(l)));
               end
               l = l + 1;
            end
         end
      end
      if (isempty(vers)),
         fprintf(2,'Version string not found!\n');
      else
         fprintf(2,'Version set to %s\n', vers);
      end
      % find the corresponding references
      ti = char(cs(ref(1)).s); ar = char(cs(ref(2)).s);
      k = 1;
      l = 1;
      ar_tmp = '';
      ti_tmp = '';
      while (k<=size(ar,1)),
          if (sum(isspace(ar(k))) == size(k,2) && ~isempty(ar_tmp)), % all space, combine
              ca(l) = cellstr([ar_tmp '-x->' ti_tmp]);
              ar_tmp = '';
              ti_tmp = '';
              l = l + 1;
          else %catenate
              ar_tmp = [ar_tmp ' ' deblank(ar(k,:))];
              ti_tmp = [ti_tmp ' ' deblank(ti(k,:))];
          end
          k = k + 1;
      end
      %Now fix up the looks of the list DO NOT CHANGE USAGE OF '-----> ' for padding, C function PrettyRef depends on them
      imx = 0;
      for (l=1:1:length(ca)),
          isep = strfind(char(ca(l)), '-x->');
          imx = max(isep, imx);
      end
      for (l=1:1:length(ca)),
          isep = strfind(char(ca(l)), '-x->');
          ca_tmp = char(ca(l));
          c1 = pad_strn(ca_tmp(1:isep), '-', imx, -1);
          c2 = ca_tmp(isep+4:length(ca_tmp));
          ca(l) = cellstr([c1 '>' c2]);
      end
      aur = char(au);
      iout = find (aur < 32 | aur > 127); %replace characters outside of basci ascii text
      aur(iout) = '-'; %dunno what to do yet....
      otr = flipud(char(ot));
      car = char(ca);
      sdecl = sprintf('char CA_EZ_REF_STR[%d][%d]', size(otr,1)+size(aur,1)+size(car,1)+10, max([size(otr,2)+15, size(aur,2)+15,size(car,2)+15]));
      sdecl2 = sprintf('char CA_EZ_VERSION_STR[%d]',length(vers)+3);
      %do someting nice
      rname  = sprintf('%s/%s', toolbox_dir, rname);
      fida = fopen(rname,'w');
      fprintf(fida, '%s = {\n',sdecl);
      fprintf(fida, '"%s",\n"%s",\n"%s",\n',otr(1,:), otr(2,:), otr(3,:));
      for (i=1:1:size(aur,1)), fprintf(fida, '"   %s",\n', aur(i,:)); end
      fprintf(fida, '"%s",\n', otr(4,:));
      for (i=1:1:size(car,1)), fprintf(fida, '"   %s",\n', car(i,:));  end
      fprintf(fida, '"%s",\n', otr(5,:));
      fprintf(fida, '" ",\n" ",\n"AFNI adaptation by",\n" Ziad S. Saad and Daniel Glen (SSCC/NIMH/NIH)",\n');
      fprintf(fida, '" Info automatically created with CA_EZ_Prep_genx.m based on se_note.m",\n');
      fclose(fida);
   end
else
   fprintf(2,'Failed to locate se_note.m\nNo new reference string created');
   return;
end

%Make NIML output files with atlas point lists
   % start with MPM (maximum probability map)
   MPMname = sprintf('%s/%s', toolbox_dir, MPMname);
   fidc = fopen (MPMname,'w');
   if (fidc < 0),
      fprintf(2,'Failed to open output NIML output file for %s\n', MPMname);
      return;
   end

% str = sprintf('# Data for atlases from Eickhoff''s SPM toolbox.\n# Automatically compiled from: %s\n# located at: %s\n by function %s\nDate: %s*/\n\n',...
%                   MPM_file(1).name, toolbox_dir, which(FuncName), date)   ;
%   fprintf(fidc,'%s', str);
%Now create MPM structure in NIML atlas point list format
   fprintf(fidc,'# -----------     MPM      ---------------------\n');
   fprintf(fidc,'# ----------- Based on: %s --------------\n', MPM_file);
   fprintf(fidc,'<atlas_point_list\n');
   fprintf(fidc,' ni_form="ni_group" >\n');
   for (i=1:1:NLbl_MPM),
      [err,PathString,FileString] = GetPath (MapMPM(i).ref, 1);

      fprintf(fidc, '<ATLAS_POINT\n');
      fprintf(fidc, '  data_type="atlas_point"\n');
      fprintf(fidc, '  STRUCT="%s"\n',deblank(MapMPM(i).name));
% Note intensity value is back to GV as it had been in the past!
%  scaleslope factor in NIFTI dataset is used and needs to be rounded off properly
      fprintf(fidc, '  VAL="%d"\n',MapMPM(i).GV);
      fprintf(fidc, '  OKEY="%d"\n',MapMPM(i).GV);
      fprintf(fidc, '  GYoAR="0"\n');
      fprintf(fidc, '  COG="0.0 0.0 0.0"\n');
      fprintf(fidc, '  />\n\n');
   end
   fprintf(fidc, '</atlas_point_list>\n');
   fclose(fidc);

   % now, repeat with the closely related PM maps (probability maps)
   PMname = sprintf('%s/%s', toolbox_dir, PMname);
   fidc = fopen (PMname,'w');
   if (fidc < 0),
      fprintf(2,'Failed to open output NIML output file for %s\n', PMname);
      return;
   end

%   str = sprintf('# Data for atlases from Eickhoff''s SPM toolbox.\n# Automatically compiled from: %s\n# located at: %s\n by function %s\nDate: %s*/\n\n',...
%                   MPM_file(1).name, toolbox_dir, which(FuncName), date);
%   fprintf(fidc,'%s', str);
%Now create PMaps structure in NIML atlas point list format
   fprintf(fidc,'# -----------     PMaps      ---------------------\n');
   fprintf(fidc,'# ----------- Based on: %s --------------\n', MPM_file);
   fprintf(fidc,'<atlas_point_list\n');
   fprintf(fidc,' ni_form="ni_group" >\n');
   for (i=1:1:NLbl_MPM),
      [err,PathString,FileString] = GetPath (MapMPM(i).ref, 1);
      fprintf(fidc, '<ATLAS_POINT\n');
      fprintf(fidc, '  data_type="atlas_point"\n');
      fprintf(fidc, '  STRUCT="%s"\n',deblank(MapMPM(i).name));
      % assume sequential structures - sub-brick 0 has STRUCT 0 with
      % SB_LABEL 0; sub-brick 1 has struct 1 with SB_LABEL 1 (0-based)
      fprintf(fidc, '  VAL="%d"\n',i-1);
      fprintf(fidc, '  OKEY="%d"\n',i-1);
      fprintf(fidc, '  GYoAR="0"\n');
      fprintf(fidc, '  COG="0.0 0.0 0.0"\n');
      fprintf(fidc, '  SB_LABEL="%s" />\n\n',deblank(RemoveExtension(FileString,'.img|.mnc|.hdr')));
   end
   fprintf(fidc, '</atlas_point_list>\n');
   fclose(fidc);



   % now, make Macrolabel NIML atlas point list
   MLname = sprintf('%s/%s', toolbox_dir, MLname);
   fidc = fopen (MLname,'w');
   if (fidc < 0),
      fprintf(2,'Failed to open output NIML output file for %s\n', MLname);
      return;
   end
   NLbl_ML = size(MapML,1);
   MaxLbl_ML = size(MapML,2)+3;

%   str = sprintf('# Data for atlases from Eickhoff''s SPM toolbox.\n# Automatically compiled from: %s\n# located at: %s\n by function %s\nDate: %s*/\n\n',...
%                   ML_file(1).name, toolbox_dir, which(FuncName), date);
%   fprintf(fidc,'%s', str);
%Now create ML (macrolabel) structure in NIML atlas point list format
   fprintf(fidc,'# -----------     Macrolabels    ---------------------\n');
   fprintf(fidc,'# ----------- Based on: %s --------------\n', ML_file(1).name);
   fprintf(fidc,'<atlas_point_list\n');
   fprintf(fidc,' ni_form="ni_group" >\n');
   for (i=1:1:NLbl_ML),
      fprintf(fidc, '<ATLAS_POINT\n');
      fprintf(fidc, '  data_type="atlas_point"\n');
      fprintf(fidc, '  STRUCT="%s"\n',deblank(MapML(i,:)));
      fprintf(fidc, '  VAL="%d"\n',i);
      fprintf(fidc, '  OKEY="%d"\n',i);
      fprintf(fidc, '  GYoAR="0"\n');
      fprintf(fidc, '  COG="0.0 0.0 0.0"\n');
      fprintf(fidc, '  />\n\n');
   end
   fprintf(fidc, '</atlas_point_list>\n\n');
   fclose(fidc);


   % now, make Left/Right Brain NIML atlas point list
   LRname = sprintf('%s/%s', toolbox_dir, LRname);
   fidc = fopen (LRname,'w');
   if (fidc < 0),
      fprintf(2,'Failed to open output NIML output file for %s\n', LRname);
      return;
   end

%   str = sprintf('# Data for atlases from Eickhoff''s SPM toolbox.\n# Automatically compiled from: %s\n# located at: %s\n by function %s\nDate: %s*/\n\n',...
%                   ML_file(1).name, toolbox_dir, which(FuncName), date);
%   fprintf(fidc,'%s', str);
%Now create LR (left/right) structure in NIML atlas point list format
   fprintf(fidc,'# -----------     LeftRight    ---------------------\n');
   fprintf(fidc,'# ----------- Based on: %s --------------\n', LRname);
   fprintf(fidc,'<atlas_point_list\n');
   fprintf(fidc,' ni_form="ni_group" >\n');
   Lst = {'Right Brain' 'Left Brain'};
   for (i=1:1:2),
      fprintf(fidc, '<ATLAS_POINT\n');
      fprintf(fidc, '  data_type="atlas_point"\n');
      fprintf(fidc, '  STRUCT="%s"\n',char(Lst(i)));
      fprintf(fidc, '  VAL="%d"\n',i);
      fprintf(fidc, '  OKEY="%d"\n',i);
      fprintf(fidc, '  GYoAR="0"\n');
      fprintf(fidc, '  COG="0.0 0.0 0.0"\n');
      fprintf(fidc, '  />\n\n');
   end
   fprintf(fidc, '</atlas_point_list>\n');
   fclose(fidc);

return;

function str = fix_string(stri)
   str = stri;
   n_str = length(str);
   %are you missing a ) ?
   i = n_str;
   broken = 0;
   closed = 0;
   while (i > 1 && ~broken),
      if (str(i) == ')'),
         closed = closed + 1;
      elseif (str(i) == '('),
         if (closed == 0),
            broken = 1;
         end
      end
      i = i - 1;
   end

   if (broken),
      i = n_str;
      fixed = 0;
      while (i > 1 && ~fixed),
         if (~isspace(str(i))),
            if (i==n_str), str = [str,')'];
            else str(i+1) = ')';
            end
            fixed = 1;
         end
         i = i - 1;
      end
   end
return;

function str = pad_with_dot(stri, ntot)
   if (nargin == 2),
      str = pad_strn(stri, ' ', ntot, -1);
   else
      str = stri;
   end
   n_str = length(str);
   i = n_str;
   while (i > 1 && isspace(str(i))),
      str(i) = '.';
      i = i - 1;
   end

return;
