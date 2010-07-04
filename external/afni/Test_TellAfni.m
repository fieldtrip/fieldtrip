%script Test_TellAfni
%
%
%
%Purpose:
%   
%   A script to demonstrate the use of the matlab AFNI driver tools (TellAfni).
%   Make sure no current AFNI session is running with the -yesplugouts option.
%
%   The script is not fancy and some steps might go by too quickly but it should 
%   be a simple read to figure it all out.
%
%Input:
%   
%   Needs the datasets distributed with AFNI's matlab library
%   http://afni.nimh.nih.gov/pub/dist/data/afni_matlab_data.tgz
%   
%Output:
%
%  Follow instructions, watch AFNI 
%   
%   
%      
%More Info :
%   
%    TellAfni
%    TellAfni_Commands
%    NewCs
%    AFNI's README.driver file
%    AFNI's plugout_drive program
%
%     Author : Ziad Saad
%     Date : Tue Dec 6 14:17:34 EST 2005
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Debug Flag
DBG = 1;

%get the directory
dirname = uigetdir(cd,'Select directory that has AFNI''s matlab demo data');
%dirname = '/Users/ziad/DownLoad/Demo_Bricks'

%check for dsets
if (exist(sprintf('%s%cARzsspgrax+orig.HEAD',dirname, filesep),'file') ~= 2),
   fprintf(2,'Error: Could not find test data in selected directory:\n%s\n', dirname);
   return;
end
%launch afni
cs = NewCs('start_afni', '', dirname);
err = TellAfni(cs);
if (err),
   fprintf(2,'Error: Failed to start AFNI in listening mode.\n');
   return;
end

%switch to relevant datsets
i = 1;
cs(i) = NewCs('Set_Anatomy', 'A', 'ARzsspgrax'); i = i + 1;
cs(i) = NewCs('open_window', '', 'axialimage', 'mont=2x2:8 keypress=v geom=500x500+800+50'); i = i+1;
err = TellAfni(cs); clear cs
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end

fprintf(1,'Sleeping for a few seconds...\n'); pause(4);
i = 1;
cs(i) = NewCs('open_window', '', 'axialimage', 'keypress=" "'); i = i+1; % stop the video with space press
cs(i) = NewCs('OPEN_PANEL', '', 'Define_Overlay'); i = i+1;
cs(i) = NewCs('Set_Function', 'A', 'ARzs_CW_avvr.DEL'); i = i + 1;
cs(i) = NewCs('See_Overlay', '', '+'); i = i + 1;
cs(i) = NewCs('SET_DICOM_XYZ', '', '-6 86 -3'); i = i+1; 
cs(i) = NewCs('SET_PBAR_SIGN', '' ,'+'); i = i + 1;
cs(i) = NewCs('SET_PBAR_NUMBER', '' ,'20'); i = i + 1;
cs(i) = NewCs('SET_SUBBRICKS', '', '-1 0 2'); i = i + 1;
cs(i) = NewCs('SET_FUNC_RANGE', '', 30); i = i + 1;
cs(i) = NewCs('SET_THRESHNEW','', 1e-9, '*p'); i = i + 1;
cs(i) = NewCs('SET_FUNC_RESAM','', 'Cu.Cu'); i = i + 1;
err = TellAfni(cs); clear cs
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end

fprintf(1,'Sleeping for a few seconds...\n'); pause(4);
i = 1;
cs(i) = NewCs('open_window', 'B', 'coronalgraph', 'geom=500x500+50+550'); i = i+1;
cs(i) = NewCs('Set_Anatomy', 'B', 'ARzs_CW_avvr+orig'); i = i+1;
cs(i) = NewCs('SET_DICOM_XYZ', 'B', '-6 86 -3'); i = i+1;
err = TellAfni(cs); clear cs
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end
 
fprintf(1,'Sleeping for a few seconds...\n'); pause(4);
i = 1;
cs(i) = NewCs('open_window', 'A', 'coronalimage', 'geom=500x500+550+750'); i = i+1; 
cs(i) = NewCs('open_window', '', 'axialimage', 'mont=1x1'); i = i+1;
err = TellAfni(cs); clear cs
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end

fprintf(1,'Sleeping for a few seconds...\n'); pause(4);
for (k=1:1:20),
   i = 2*k-1;
   cs(i) = NewCs('PBAR_ROTATE', '', '+'); i = i+1;
   fnm = sprintf('Rot_%s.jpg',pad_strn(sprintf('%d',k), '0', 2, 1));
   unix(sprintf('rm %s', fnm));
   cs(i) = NewCs('SAVE_JPEG', '', 'coronalimage', fnm); 
end
err = TellAfni(cs); clear cs
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end

%load then show the images written to disk
for (i=1:1:20),
   fnm = sprintf('Rot_%s.jpg',pad_strn(sprintf('%d',i), '0', 2, 1));
   ts(i).im = imread(fnm);
end
figure(1); clf;
for (i=1:1:200),
   imshow(ts(rem(i,20)+1).im); drawnow
end   


input ('All done, hit "enter" to quit\n','s');
err = TellAfni(NewCs('Quit'));
if (err),
   fprintf(2,'Error: Failed telling AFNI.\n');
   return;
end
