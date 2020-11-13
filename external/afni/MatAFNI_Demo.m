%script MatAFNI_Demo
%
%
%
%Purpose:
%   Demonstrate the use of the matlab AFNI library
%   The script is echoed onto the screen but it is best
%   you read through it with a text editor.
%   
%   
%Input:
%   You must have the data sets 
%   ARzs_CW_avvr+orig, ARzs_CW_avvr.DEL+orig and ARzsspgrax+orig
%   in the working directory
%
%Output:
%   
%
%   
%   
%      
%Key Terms:
%   
%More Info :
%   
%   
%   
%
%     Author : Ziad Saad
%     Date : Tue Aug 28 16:21:54 PDT 2001
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland

echo on

clear all

%check for data sets before proceeding 
	if (exist('ARzs_CW_avvr+orig.HEAD') ~= 2 | exist('ARzs_CW_avvr+orig.BRIK') ~= 2),
		fprintf ('\nERROR: Missing ARzs_CW_avvr+orig data.\nCtrl+c to Abort.\n');
		while (1), 
			pause;
		end
	end

	if (exist('ARzs_CW_avvr.DEL+orig.HEAD') ~= 2 | exist('ARzs_CW_avvr.DEL+orig.BRIK') ~= 2),
		fprintf ('\nERROR: Missing ARzs_CW_avvr.DEL+orig data.\nCtrl+c to Abort.\n');
		while (1), 
			pause;
		end
	end
	if (exist('ARzsspgrax+orig.HEAD') ~= 2 | exist('ARzsspgrax+orig.BRIK') ~= 2),
		fprintf ('\nERROR: Missing ARzsspgrax+orig data.\nCtrl+c to Abort.\n');
		while (1), 
			pause;
		end
	end

%%%%%%%%%%%%%% Read time series and access time series %%%%%%%%%%%%%%%%%%%%%%%%
	fprintf (1,'\tpatience...\n'); 
	%read in a brick and display voxel with AFNI i j k: 29, 33, 3 
		BrikName = 'ARzs_CW_avvr+orig';
		Vox_Ind = [29 33 3];

	%Do it using matrix format for V.
		Opt.Format = 'matrix';
		[err, Vm, Info, ErrMessage] = BrikLoad (BrikName, Opt); 
	%plot it
 		figure(1); clf; subplot (211);
		plot (squeeze(Vm(Vox_Ind(1)+1,Vox_Ind(2)+1, Vox_Ind(3)+1, :)));


	%Do it using vector format for V. Note that this example does not highlight
	%the advantage of using the vector format. 
		Opt.Format = 'vector';
		[err, Vv, Info, ErrMessage] = BrikLoad (BrikName, Opt); 
	%figure out the AFNI 1D index of the voxel
		[err, Indx] = AfniXYZ2AfniIndex (Vox_Ind, Info.DATASET_DIMENSIONS(1), Info.DATASET_DIMENSIONS(2));
	%plot it. Remember that indexing in matlab is always augmented by 1 relative to AFNI's 
 		figure(1); subplot (212);
		plot (Vv(Indx+1,:));
		plotsign2(1);

	drawnow
	msgbox(sprintf('Showing Time Series From %s.\nFrom AFNI''s viewer, use jump to:  %g %g %g to see the same time series',...
	  BrikName, Vox_Ind(1),Vox_Ind(2), Vox_Ind(3))); 
	fprintf (1,'\nHit Enter to Proceed with Demo  (ctrl+c to abort)...\n'); pause
	
%%%%%%%%%%%%%% Read Anatomy and grab a slice %%%%%%%%%%%%%%%%%%%%%%%%
	fprintf (1,'\tpatience...\n'); 
	%read in anatomy brick
		BrikNameAnat = 'ARzsspgrax+orig';
		[err, Vanat, InfoAnat, ErrMessage] = BrikLoad (BrikNameAnat);	
	
	%Grab a few slices in the same plane
		OptDisp.plane = 'Ax';
		OptDisp.iSlc = [107 79 45]; %Those are AFNI indices, see help GetAfniSlice for more info
		[err, slc, slcDisp] = GetAfniSlice (Vanat, InfoAnat, OptDisp);
		figure(2); clf
		subplot 321; imagesc (slc.M(:,:,1)); title (sprintf('Ax #%g, As In Brick', OptDisp.iSlc(1)));
		subplot 322; imagesc (slcDisp.M(:,:,1)); title (sprintf('Ax #%g, As Displayed', OptDisp.iSlc(1)));
		subplot 323; imagesc (slc.M(:,:,2));title (sprintf('Ax #%g, As In Brick', OptDisp.iSlc(2)));
		subplot 324; imagesc (slcDisp.M(:,:,2));	title (sprintf('Ax #%g, As Displayed', OptDisp.iSlc(2)));
		subplot 325; imagesc (slc.M(:,:,3));title (sprintf('Ax #%g, As In Brick', OptDisp.iSlc(3)));
		subplot 326; imagesc (slcDisp.M(:,:,3));	title (sprintf('Ax #%g, As Displayed', OptDisp.iSlc(3)));
		colormap gray %you need to use some scaling to get the display to work best
		plotsign2(2);
		
		drawnow
		msgbox(sprintf('Showing Anatomy Brick %s.\nFrom AFNI''s axial viewer, use slider to view slices %g, %g and %g that are displayed in the right column of the figure.',...
		 BrikNameAnat, OptDisp.iSlc(1),OptDisp.iSlc(2), OptDisp.iSlc(3))); 
		fprintf (1,'\nHit Enter to Proceed with Demo  (ctrl+c to abort)...\n'); pause
		
%%%%%%%%%%%%%% Read Functional brick and grab a slice triplet %%%%%%%%%%%%%%%%%%%%%%%%
	fprintf (1,'\tpatience...\n'); 
	%grab a triplet of slices from .DEL functional Brick
		BrikNameFunc = 'ARzs_CW_avvr.DEL+orig';
		[err, Vfunc, Infofunc, ErrMessage] = BrikLoad (BrikNameFunc);	
		
		DM = AFNI_SliceDispManip (Infofunc);
		OptDispFunc.iSlc = [31 32 6];
		OptDispFunc.index = 1;
		[err, slc] = GetAfniSliceTriplet (Vfunc, Infofunc, DM, OptDispFunc);
		figure(3); clf
		subplot 321; imagesc (slc(1).M(:,:)); title (sprintf('Ax #%g, As In Brick', OptDispFunc.iSlc(1))); 
		subplot 322; imagesc (slc(1).Mdisp(:,:)); title (sprintf('Ax #%g, As Displayed', OptDispFunc.iSlc(1)));
		subplot 323; imagesc (slc(2).M(:,:));title (sprintf('Sa #%g, As In Brick', OptDispFunc.iSlc(2)));
		subplot 324; imagesc (slc(2).Mdisp(:,:));	title (sprintf('Sa #%g, As Displayed', OptDispFunc.iSlc(2)));
		subplot 325; imagesc (slc(3).M(:,:));title (sprintf('Co #%g, As In Brick', OptDispFunc.iSlc(3)));
		subplot 326; imagesc (slc(3).Mdisp(:,:));	title (sprintf('Co #%g, As Displayed', OptDispFunc.iSlc(3)));
		colormap gray %you need to use some scaling to get the display to work best
		plotsign2(3);
		
		drawnow
		msgbox(sprintf('Showing Slice Triplet from Functional Brick %s.\nFrom AFNI''s axial viewer, use slider to view slices %g, %g and %g that are displayed in the right column of the figure.\nSet Func underlay and set Func sub-brick to #%g',...
		 BrikNameFunc, OptDispFunc.iSlc(1),OptDispFunc.iSlc(2), OptDispFunc.iSlc(3), OptDispFunc.index)); 
		fprintf (1,'\nHit Enter to Proceed with Demo (ctrl+c to abort)...\n'); pause


%%%%%%%%%%%%%%%%% Show a cool graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		fprintf (1,'\tpatience...\n'); 
		%Extract the delay and cross correlation coefficient sub-bricks
		%see Infofunc.BRICK_LABS for indexing info
		Vdel = squeeze(Vfunc(:,:,:,1)); %that's not necessary but it makes annotation and referencing clear
		Vxcor = squeeze(Vfunc(:,:,:,3));		
		
		%find the cross correlation coefficients < 0.5
		inoact = find(Vxcor < 0.5);
		
		%mask those voxels with  -1
		Vdel(inoact) = -1;
		
		%make sure no voxels having a delay outside [0 .. 40] are masked
		iout = find(Vdel < 0 | Vdel > 40);
		Vdel(iout) = -1; %Note that for efficiency, the previous two steps can be combined into 1
		
		%also mask the Cross Correlation Subbrick for later use
		Vxcor(iout) = 0;
		
		%Extract axial slice 29
		OptFunc2.iSlc = 29;
		OptFunc2.plane = 'Ax';
		[err, slc, slcDisp] = GetAfniSlice (Vdel, Infofunc, OptFunc2);
		
		%interpolate for pretty pictures
		SlcShow = interp2(slcDisp.M, 3, 'cubic');
		
		%display a surface mesh of delay values for slice 
		figure(4); clf
		subplot (211); 
		hs = surf(SlcShow, 'EdgeColor', 'none'); title ('Surface Plot of Delay Values');
		subplot (223);
		hc = contour (SlcShow); title ('Contour Plot');
		subplot (224);
		ipos = find(Vdel >= 0);
		hist(Vdel(ipos),20); title ('Histogram');
		plotsign2(4);
		drawnow;
		fprintf (1,'\nHit Enter to Proceed with Demo (ctrl+c to abort)...\n'); pause

%%%%%%%%%%%%%%%%% Writing Out Functional Data %%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Write out delay brick with only the masked voxels showing
	
	fprintf (1,'\tpatience...\n'); 
	%We need to setup the header info for the new data set to write out
	%Copy the header info from the .DEL brick
	InfoDelOut = Infofunc;
	%modify the necessary fields 
	%YOU MUST READ THE HELP for the function WriteBrik and BrikInfo
	InfoDelOut.RootName = ''; %that'll get set by WriteBrik
	InfoDelOut.DATASET_RANK(2) = 2; %two sub-bricks
	InfoDelOut.BRICK_TYPES = [1 1]; %store data as shorts
	InfoDelOut.BRICK_STATS = []; %automatically set
	InfoDelOut.BRICK_FLOAT_FACS = []; %automatically set
	InfoDelOut.BRICK_LABS = 'Masked Delay~Masked Cross Correlation';
	InfoDelOut.IDCODE_STRING = ''; %automatically set
	InfoDelOut.BRICK_STATAUX = [1 2 3 160 2 2];
	
	OptDelOut.Scale = 1;
	OptDelOut.Prefix = 'Demo1_func';
	OptDelOut.verbose = 0;
	!rm Demo1_func*
	%To write the two sub-bricks, I placed them in an Nx2 matrix [Vdel(:) Vxcor(:)]
	[err, ErrMessage, InfoDelOut] = WriteBrik ([Vdel(:) Vxcor(:)], InfoDelOut, OptDelOut);

	if (err),
		errordlg(sprintf('%s\nAbort now (ctrl+c).', ErrMessage));
		while (1),
			pause;
		end
	end
	
	msgbox(sprintf('Data set %s+orig was written to disk, check it out with AFNI.',...
		 OptDelOut.Prefix)); 

	fprintf (1,'\nHit Enter to Proceed with Demo (ctrl+c to abort)...\n'); pause

%%%%%%%%%%%%%%%%% Writing Out time series %%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf (1,'\tpatience...\n'); 
	%Write out random time series where functional voxels passed the threshold above
	%find voxel indices that are not masked
	[i_in] = setdiff([1:1:prod(Info.DATASET_DIMENSIONS(1:3))],iout);
	%create that many random time series  
	Mrand = randn(length(i_in), Info.DATASET_RANK(2)); 
	%create an empty volume
	Vnew = zeros(prod(Info.DATASET_DIMENSIONS(1:3)), Info.DATASET_RANK(2));
	%replace zeros with random time series
	Vnew(i_in,:) = Mrand;
	
	%prepare Header for writing
	InfoNewTSOut = Info;
	InfoNewTSOut.RootName = '';
	InfoNewTSOut.BRICK_STATS = []; %automatically set
	InfoNewTSOut.BRICK_FLOAT_FACS = []; %automatically set
	InfoNewTSOut.IDCODE_STRING = ''; %automatically set
	 			
	OptTSOut.Scale = 1;
	OptTSOut.Prefix = 'Demo1_TS';
	OptTSOut.verbose = 1;
	!rm Demo1_TS*
	%write it
	[err, ErrMessage, InfoNewTSOut] = WriteBrik (Vnew, InfoNewTSOut, OptTSOut);
	msgbox(sprintf('Data set %s+orig was written to disk, check it out with AFNI.',...
		 OptTSOut.Prefix)); 
	
	if (err),
		errordlg(sprintf('%s\nAbort now (ctrl+c).', ErrMessage));
		while (1),
			pause;
		end
	end
	
	fprintf (1,'\nDemo Done\n');

echo off				
