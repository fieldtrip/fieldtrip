function [Cs] = rgbdectohex (Mrgb,strg)
%
%  Mode 1:
%        [Cs] = rgbdectohex (Mrgb,[string])
%  Mode 2:
%        rgbdectohex 
%
% Mode 1:
%This function takes as an input the rgb matrix Mrgb (size is Nx3)
% where each row represents the rgb gun values of a color.
% gun values should be integers between 0 and 255
%
%the optional 'string' is placed before the #hex representation 
%in a way VERY suitable for afni's .Xdefaults file stuff ...
%
% remember the magic command : xrdb -merge <colors filename>
% 
% the function cyclmat helps you fix your color maps, check it out...
%
%The hex values are all padded to two characters, looks nice.
%
% The function first displays the map, then  gives you the option of 
%   writing out the results:
%     to an ascii file that can be (if you specified the right strg parameter) used 
%        directly in .Xdefaults
%     to an ascii file that contains RGB values
%     to an ascii file containing the definition of the colormap in a C-
%        syntax that can be added to pbar.c and added with PBAR_define_bigmap( char *mapcmd );
% 
%the result is written to stdout in a format used 
%by .Xdefault files, to make importing them to
%.Xdefaults easy, just cut and paste.
%
%example: >> load hues_rygbr20 (this file is an ascii list of 3 integers per line specifying rgb colours)
%         >> Mcyc = cyclmat (hues_rygbr20,-1,13); (changes the order of the loaded color file..)
%         >> rgbdectohex (Mcyc,'AFNI*ovdef'); this displays the colour map, and asks if you want the results written out to a file
%         
% If you sqaved the results to a file called junkmap
% from command line do : xrdb -merge junkmap and the colormap is set for afni to read.
%
% If you use the return parameters [Cs, sall] then you'll get a command structure vector
%  that tells AFNI (TellAfni(Cs)) to load the newly created colorscale and switch to it.
%  The colorscale is named string if one is supplied.  
%
% Mode 2:
% Interactive mode for showing RGB colors. Enter RGB values to see color.
%
% see also MakeColorMap, TellAfni, ROIcmap
%
%		Ziad Saad Nov 26 97/  Dec 4 97/ Jan 06

if (nargout), Cs = []; sall = ''; end

if (nargin == 2),
	son = 1;
else
	son = 0;
   strg = '';
end

if (nargin == 0),
   while (1),
      strgb = input ('Enter rgb values (nothing to exit):','s');
      if (isempty(strgb)) return; end
      Mrgb = str2num(strgb);
      ShowRGBcol(Mrgb);
      fprintf(1,'        Color: %s\n', RGBtoXhex(Mrgb, 0, ''));
   end
end

if (size (Mrgb,2) ~= 3)
	fprintf(2,'rgbdectohex : Wrong Mrgb matrix size')
	Mhex = -1;
	return;
end

if (max(abs(round(Mrgb(:))-Mrgb(:))) > 0.0001),
   fprintf(2,'rgbdectohex : RGB values are not integers.\n');
   Mhex = -1;
	return;
end

if (max(Mrgb(:)) < 1.1),
   fprintf(2,'rgbdectohex : Maximum RGB value < 1.1. Suspecting RGB to be between 0 and 1\nBe sure that color values are integers that can range between 0 and 255.')
	input ('Hit enter to continue','s');
end

if (min(Mrgb(:)) < 0 | max(Mrgb(:)) > 255),
   fprintf(2,'rgbdectohex : RGB value must range between 0 and 255.');
   return;
end

ShowRGBcol(Mrgb);

for (i=1:1:size(Mrgb,1)),
		fprintf(1,'%s\n',RGBtoXhex(Mrgb(i,:), son, strg));
end %i

chc = input ('Wanna write this to disk ? (y/n)','s');
if (chc == 'y'),
	chc2 = input ('Write rgb version of map too ? (y/n)','s');
	chc3 = input ('Write C-friendly version of map to include in AFNI''c code ? (y/n)','s');
	rep = 1;
   strout = strg;
	while (rep == 1),
		%strout = input ('Enter filename :','s');
		rep = filexist (strout);
		if (rep == 1),
			fprintf (2,'rgdtodec : file %s exists, enter another name\n\a', strout);
		end
	end
	strout2 = sprintf ('%s.rgb',strout);
	strout3 = sprintf ('%s.love.c',strout);
	fid = fopen (strout,'wt');
	if (chc2 == 'y'),
		fid2 = fopen (strout2,'wt');
	end
   if (chc3 == 'y'),
		fid3 = fopen (strout3,'wt');
	end
	
	if (fid == -1 | (chc2 == 'y' & fid2 == -1) | (chc3 == 'y' & fid3 == -1)),
		fprintf (2,'rgdtodec : Could not open %s or %s  %s file for write operation\n',strout,strout2, strout3);
	end
	
   if (chc3 == 'y'),
         fprintf (fid3,'static char %s_CMD[] = {  \n      "%s "\n      "',  upper(strg), strg);
   end
	for (i=1:1:size(Mrgb,1)),
		s1 = pad_strn (lower(dec2hex(Mrgb(i,1))),'0',2,1);
		s2 = pad_strn (lower(dec2hex(Mrgb(i,2))),'0',2,1);
		s3 = pad_strn (lower(dec2hex(Mrgb(i,3))),'0',2,1);
	
		if (son == 1),
			cst = sprintf ('%g',i);
			cpd = pad_strn (cst,'0',2,1);
			fprintf (fid,'%s%s:\t#%s%s%s\n',strg,cpd,s1,s2,s3);	
		else
			fprintf (fid,'#%s%s%s\n',s1,s2,s3);
		end
		if (chc2 == 'y'),
			fprintf (fid2,'%g %g %g\n',Mrgb(i,1),Mrgb(i,2),Mrgb(i,3));
		end
      if (chc3 == 'y'),
         if (rem(i, 5)==1 & i > 1),
            fprintf (fid3,'"\n      "');
         end
         fprintf (fid3,'#%s%s%s ',s1,s2,s3);
      end
	end %i
	fclose (fid);
	if (chc2 == 'y'),
		fclose (fid2);
	end
   if (chc3 == 'y'),
      fprintf (fid3,'"\n};');
		fclose (fid3);
	end
end

if (nargout),
   sall = '';
   for (i=1:1:size(Mrgb,1)),
      s1 = pad_strn (lower(dec2hex(Mrgb(i,1))),'0',2,1);
		s2 = pad_strn (lower(dec2hex(Mrgb(i,2))),'0',2,1);
		s3 = pad_strn (lower(dec2hex(Mrgb(i,3))),'0',2,1);
      
      sall = sprintf('%s#%s%s%s ',sall, s1,s2,s3); 
   end
   if (son == 0), strg = 'rgbdectohex'; end
   Cs = NewCs('DEFINE_COLORSCALE', '', strg, sall);
   Cs(2) = NewCs('SET_PBAR_ALL', '' , '+99', sprintf('1.0 %s', strg));
end

return;

function ShowRGBcol(Mrgb)
figure (1);
colormap (Mrgb./255);
subplot 211;
image ([1:1:length(Mrgb(:,1))]);
title(sprintf('X11 Color in Hex: %s', RGBtoXhex(Mrgb, 0, '')));

%subplot 212;
%pie (ones(1,length(Mrgb(:,1))));

return;

function sret = RGBtoXhex(rgb, son, strg),
      s1 = pad_strn (lower(dec2hex(rgb(1))),'0',2,1);
		s2 = pad_strn (lower(dec2hex(rgb(2))),'0',2,1);
		s3 = pad_strn (lower(dec2hex(rgb(3))),'0',2,1);
	
		if (son == 1),
			cst = sprintf ('%g',i);
			cpd = pad_strn (cst,'0',2,1);
			sret = sprintf ('%s%s:\t#%s%s%s',strg,cpd,s1,s2,s3);	
		else
			sret = sprintf ('#%s%s%s',s1,s2,s3);
		end
return;
