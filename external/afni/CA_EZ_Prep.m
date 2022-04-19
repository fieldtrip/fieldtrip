function [MapMPM, MapML, aur] = CA_EZ_Prep()
% A function to process a new Eickhoff, Amunts and Zilles CA toolbox
% The function creates new versions of AFNI's source files
% thd_ttatlas_CA_EZ.c and thd_ttatlas_CA_EZ.h .
% The new files created are called thd_ttatlas_CA_EZ-auto.c and thd_ttatlas_CA_EZ-auto.h
% must be inspected then moved and renamed into afni's src directory.
%
% See SUMA/Readme_Modify.log for info on sequence of execution
% search for: + How you install a new Zilles, Amunts, Eickhoff SPM toolbox:%
%
% See also scripts:
%  @Prep_New_CA_EZ
%  @Compare_CA_EZ
%  @Create_suma_tlrc.tgz
%  @Create_ca_ez_tlrc.tgz
%  @DistArchives
% ZSS SSCC Feb 06

FuncName = 'CA_EZ_Prep';
MapMPM = [];
MapML = [];

toolbox_dir = '/Volumes/afni/home4/users/ziad/Programs/matlab/spm2/toolbox/Anatomy';

if (exist(toolbox_dir) ~= 7),
   fprintf(2,'Toolbox directory %s not found\nPick a new one:', toolbox_dir);
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
   cd (curdir);
   [err,pt] = GetPath(toolbox_dir);
   [err,sout] = unix(sprintf('ls -l %s', pt));
   c = input(sprintf('Found toolbox here: %s\nDirectory Listing:\n%s\nEnter "y" to use it, anything else to quit.\n', toolbox_dir, sout),'s');
   if (isempty(c) || ( c(1) ~= 'y' && c(1) ~= 'Y' ) ),
      return;
   end
end

   fprintf(1,'Using toolbox directoy %s...\n', toolbox_dir);

%First get the MPM info
   prf = sprintf('%s%cAllAreas*MPM.mat', toolbox_dir, filesep);
   [err, ErrMessage, MPM_file] = zglobb ({prf});
   if (size(MPM_file,1) ~= 1),
      fprintf(2,'Could not find unique MPM map list\n', toolbox_dir);
      for (i=1:1:size(MPM_file,1)), MPM_file(i), end
      return;
   end

   %load the MPM map structure
   MapMPM = load(MPM_file(1).name);
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
      fprintf(2,'Could not find unique ML map list\n', toolbox_dir);
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
   cname = 'thd_ttatlas_CA_EZ-auto.c';
   hname = 'thd_ttatlas_CA_EZ-auto.h';
   rname = 'thd_ttatlas_CA_EZ-ref.h';

%Do the references
% a horrible mess of an if block below. One hopes for a better solution someday
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
         tmp = get(h(i),'String')
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
         return;
      else
         fprintf(2,'Version set to %s\n', vers);
      end
      % find the corresponding references
      ti = char(cs(ref(1)).s); ar = char(cs(ref(2)).s);
      if (nref ~= 2 || (size(ar,1) ~= size(ti,1))),
         fprintf(2,'Warning:\nUnexpected number of ref strings or some mismatch (%d, %d, %d)\nYou have to edit thd_ttat', nref, size(ar,1), size(ti,1));
         sdecl = sprintf('char CA_EZ_REF_STR[128][256]');
         sdecl2 = sprintf('char CA_EZ_VERSION_STR[128]');
         %return;
      else
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
         fida = fopen(rname,'w');
         fprintf(fida, '%s = {\n',sdecl);
         fprintf(fida, '"%s",\n"%s",\n"%s",\n',otr(1,:), otr(2,:), otr(3,:));
         for (i=1:1:size(aur,1)), fprintf(fida, '"   %s",\n', aur(i,:)); end
         fprintf(fida, '"%s",\n', otr(4,:));
         for (i=1:1:size(car,1)), fprintf(fida, '"   %s",\n', car(i,:));  end
         fprintf(fida, '"%s",\n', otr(5,:));
         fprintf(fida, '" ",\n" ",\n"AFNI adaptation by",\n" Ziad S. Saad (saadz@mail.nih.gov, SSCC/NIMH/NIH)",\n');
         fprintf(fida, '" Info automatically created with CA_EZ_Prep.m based on se_note.m",\n');
         fprintf(fida, '""};/* Must be the only empty string in the array*/\n'); %Must be the only empty string in the array
         fprintf(fida, '%s = { "%s" };\n', sdecl2, vers);
         fclose(fida);
      end
   end
else
   fprintf(2,'Failed to locate se_note.m\nNo new reference string created');
   return;
end

%Prep C output files

   fidc = fopen (cname,'w');
   if (fidc < 0),
      fprintf(2,'Failed to open output .c file\n');
      return;
   end
   fidh = fopen (hname,'w');
   if (fidh < 0),
      fprintf(2,'Failed to open output .h file\n');
      return;
   end

   str = sprintf('/*! Data for atlases from Eickhoff''s SPM toolbox.\nAutomatically compiled from: %s\n located at: %s\n by function %s\nDate: %s*/\n\n',...
                   MPM_file(1).name, toolbox_dir, which(FuncName), date);
   fprintf(fidh,'%s', str);
   fprintf(fidc,'%s', str);

%Add the references string file
   fprintf(fidc,'/*! Leave the reference string in a separate file\nfor easy script parsing.*/\n#include "%s"\n\n', rname);

%Now do specifics to .h file
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

   fprintf(fidh,'/* ----------- Macro Labels --------------------- */\n');
   fprintf(fidh,'/* ----------- Based on: %s -------------*/\n', ML_file(1).name);
   fprintf(fidh,'#define ML_EZ_COUNT   %d\n\n', NLbl_ML);
   fprintf(fidh,'extern ATLAS_point ML_EZ_list[ML_EZ_COUNT] ;\nextern char * ML_EZ_labels[ML_EZ_COUNT] ;\n');
   fprintf(fidh,'extern int ML_EZ_labeled ;\nextern int ML_EZ_current ;\n');
   fprintf(fidh,'/* ----------- Left Right   --------------------- */\n');
   fprintf(fidh,'/* ---- Based on my understanding -------------- */\n');
   fprintf(fidh,'#define LR_EZ_COUNT   3\n\n');
   fprintf(fidh,'extern ATLAS_point LR_EZ_list[LR_EZ_COUNT] ;\nextern char * LR_EZ_labels[LR_EZ_COUNT] ;\n');
   fprintf(fidh,'extern int LR_EZ_labeled ;\nextern int LR_EZ_current ;\n\n');
   fprintf(fidh,'/* -----------     MPM      --------------------- */\n');
   fprintf(fidh,'/* ----------- Based on: %s --------------*/\n', MPM_file(1).name);
   fprintf(fidh,'#define CA_EZ_COUNT   %d\n', NLbl_MPM);
   fprintf(fidh,'#define CA_EZ_MPM_MIN 100  /*!< minimum meaningful value in MPM atlas */\n');
   fprintf(fidh,'extern ATLAS_point CA_EZ_list[CA_EZ_COUNT] ;\nextern char * CA_EZ_labels[CA_EZ_COUNT] ;\n');
   fprintf(fidh,'extern int CA_EZ_labeled ;\nextern int CA_EZ_current ;\n\n');
   fprintf(fidh,'/* -----------     Refs      --------------------- */\n');
   fprintf(fidh,'/* ----------- Based on se_note.m --------------*/\n');
   fprintf(fidh,'extern %s;\n', sdecl);
   fprintf(fidh,'extern %s;\n', sdecl2);


%first create ML structure
   fprintf(fidc,'/* ----------- Macro Labels --------------------- */\n');
   fprintf(fidc,'/* ----------- Based on: %s -------------*/\n', ML_file(1).name);
   fprintf(fidc,'ATLAS_point ML_EZ_list[ML_EZ_COUNT] = {\n');
   for (i=1:1:size(MapML,1)),
      %pad string by dots
      fprintf(fidc,'   { %-3d , "%s", 0, 0, 0, 0, ""}', i, pad_with_dot(MapML(i,:), 50));
      if (i<size(MapML,1)) fprintf(fidc,',\n'); else fprintf(fidc,'\n'); end
   end
   fprintf(fidc,'};\n\n');

%Now create MPM structure
   fprintf(fidc,'/* -----------     MPM      --------------------- */\n');
   fprintf(fidc,'/* ----------- Based on: %s --------------*/\n', MPM_file(1).name);
   fprintf(fidc,'ATLAS_point CA_EZ_list[CA_EZ_COUNT] = { \n');
   for (i=1:1:NLbl_MPM),
      [err,PathString,FileString] = GetPath (MapMPM(i).ref, 1);
      fprintf(fidc,'   { %-3d, "%s", 0, 0, 0, 0, "%s" }', ...
                        MapMPM(i).GV, pad_with_dot(MapMPM(i).name,40), pad_with_dot(RemoveExtension(FileString,'.img|.mnc|.hdr'), 27));
      if (i<NLbl_MPM) fprintf(fidc,',\n'); else fprintf(fidc,'\n'); end
   end
   fprintf(fidc,'};\n\n');

%Now create LR structure
   fprintf(fidc,'/* ----------- Left Right   --------------------- */\n');
   fprintf(fidc,'/* ---- Based on my understanding -------------- */\n');
   fprintf(fidc,'ATLAS_point LR_EZ_list[LR_EZ_COUNT] = {\n');
   Lst = ['Non-Brain...'; 'Right Brain.'; 'Left Brain..'];
   for (i=1:1:3),
      fprintf(fidc,'   { %-3d, "%s", 0, 0, 0, 0, "" }', ...
                        i-1, Lst(i, :));
      if (i<3) fprintf(fidc,',\n'); else fprintf(fidc,'\n'); end
   end
   fprintf(fidc,'};\n\n');

fclose(fidh); fclose(fidc);


if (exist(rname) == 2),
   lst = sprintf('%s, %s and %s', cname, hname, rname);
else
   lst = sprintf('%s and %s (%s was not created!)', cname, hname, rname);
end

fprintf(1,'\nThe files:\n %s\n in %s\nare meant to replace\n thd_ttatlas_CA_EZ.c, thd_ttatlas_CA_EZ.h and %s\nin AFNI''s src code.\n',...
            lst, pwd, rname) ;


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
