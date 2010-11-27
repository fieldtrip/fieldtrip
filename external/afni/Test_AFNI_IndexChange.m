%script Test_AFNI_IndexChange
%
%
%
%Purpose:
%   
%   
%   
%Input:
%   
%   
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
%     Date : Fri Sep 8 14:02:23 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Debug Flag
DBG = 1;

clear all

%create a colormap
OptMap.Range = 1; OptMap.SkipLast = 0; OptMap.Showme = 0; 
[err, ColMap] = MakeColorMap ([0 0 0; 1 1 1], 256, OptMap);
colormap(ColMap); 

%load the brik
[err, V, Info] = BrikLoad('ARzs_CW_avvr+orig.HEAD'); 

Opt.plane = 'Ax';
Opt.iSlc = [2 11]; %as shown in the AFNI slider bar
[err, slc, slcDisp] = GetAfniSlice (V, Info, Opt);

OptDisp.handle = 1;
OptDisp.subplot = [2 1];
OptDisp.colrange = [2 98];
OptDisp.plane = Opt.plane;
OptDisp.Info = Info;
OptDisp.iSlc = Opt.iSlc;


[err] = DisplayAfniSlice (slcDisp, OptDisp);

	fprintf(1,'\nLeft click at 4 different locations in the top image.\n');
	[I, J] = ginput (4); %choose 4 points on slice 1
	Iorig = [round(I(:)-1) round(J(:)-1) OptDisp.iSlc(1).*ones(size(I(:)))]


	%show these points
	slcDisp.M(Iorig(1,2),Iorig(1,1) ,1) = 256;
	slcDisp.M(Iorig(2,2),Iorig(2,1),1) = 256;
	slcDisp.M(Iorig(3,2),Iorig(3,1),1) = 256;
	slcDisp.M(Iorig(4,2),Iorig(4,1),1) = 256;
	
[err] = DisplayAfniSlice (slcDisp, OptDisp);

	Iorig
  [err, Itrans] = AFNI_IndexChange (Info, Iorig, 'D2A') ;
  Itrans
  [err, Itrans] = AFNI_IndexChange (Info, Itrans, 'A2D') ;
   Itrans
                                                                                                                          
