function test_bug1826

% at this moment the test script does not yet work, but we don't want the automatic regression testing to flag it as failure
%return;
% Uses the linear tranformation algorithms of FSL to register a T2 image
% and DTI image to T1 space. It also rotates the gradient direction vectors
% .bvec with the rotational part of the transformation matrix.
% 4 different types of interpolation options:
% {trilinear,nearestneighbour,sinc,spline}
% and 5 different types of cost functions can be used:
%{mutualinfo,corratio,normcorr,normmi,leastsq}
%please do not forget to modify the fsl version as in your system

% TEST test_bug1826
% TEST ft_read_mri ft_write_mri ft_volumerealign

%%
%clear all;
subjectT1  = '/home/common/matlab/fieldtrip/data/test/bug1826/T1.nii.gz';%'~/meeting_simbio/test_bug1826_data/test_bug1826_data/T1.nii.gz';
subjectT2  = '/home/common/matlab/fieldtrip/data/test/bug1826/T2.nii.gz';%'~/meeting_simbio/test_bug1826_data/test_bug1826_data/T2.nii.gz';
subjectDTI = '/home/common/matlab/fieldtrip/data/test/bug1826/DTI.nii';%'~/meeting_simbio/test_bug1826_data/test_bug1826_data/DTI.nii';

% load in the data
T1  = ft_read_mri(subjectT1);
T2  = ft_read_mri(subjectT2);
DTI = ft_read_mri(subjectDTI);

figure; ft_plot_ortho(T1.anatomy); title('T1 before aligning')
figure; ft_plot_ortho(T2.anatomy); title('T2 before aligning')
figure; ft_plot_ortho(DTI.anatomy(:,:,:,1));title('DTI before aligning')

% they are now in memory, which would be the normal `starting point in a
% fieldtrip analysis pipeline

% INSERT THE NEW CODE HERE, THIS PROBABLY INVOLVES
% step 1) write them to disk, perhaps using ft_write_mri
% step 2) call external code, e.g. FSL (using a system call) or SPM
% step 3) read the results from disk

% define the cfg for ft_volumerealign
cfg          = [];
cfg.method   = 'volume';
cfg.fsl.path = '/opt/fsl/bin'; % '/opt/fsl_5.0/bin'; % fsl_5.0 only works on high mentats due to libraries
cfg.fsl.dof  = 6; % rigid body

interpmethods = {'nearestneighbours', 'sinc', 'trilinear'};
costfuns      = {'mutualinfo', 'corratio', 'normcorr', 'normmi', 'leastsq'};

% FIXME: not all combinations of interpmethods and cost functions make
% sense! Between modality considerations etc

%T2->T1
for k = 1:numel(interpmethods)
  for m = 1:numel(costfuns)
    cfg.fsl.interpmethod = interpmethods{k};
    cfg.fsl.costfun      = costfuns{m};
    T2aligned{k,m}       = ft_volumerealign(cfg, T2, T1);
  end
end

cfg.fsl.dof = 12;

%DTI->T1
for k = 1:numel(interpmethods)
  for m = 1:numel(costfuns)
    cfg.fsl.interpmethod = interpmethods{k};
    cfg.fsl.costfun      = costfuns{m};
    DTIaligned{k,m}      = ft_volumerealign(cfg, DTI, T1);
  end
end



%% FSL parameters

%cfg.fsl.version='/opt/fsl_4.1/bin/';% write the fsl version that you are using. The version I use is 4.1

cfg.fsl.path='/opt/fsl_4.1/bin';
%cfg.fsl.path='/opt/fsl_4.1/bin'; % should be tested if it works with different versions of FSL
cfg.fsl.T2cost='mutualinfo';% mutualinfo,corratio,normcorr,normmi,leastsq
cfg.fsl.T2dof='6';
%cfg.fsl.T2interp='spline';%trilinear,nearestneighbour,sinc,spline
%cfg.fsl.T2interp='trilinear';
cfg.fsl.T2interp='sinc';

cfg.fsl.DTIcost='mutualinfo';
cfg.fsl.DTIdof='6';
cfg.fsl.DTIinterp='spline';


%% Use until dot for naming

dots = strfind(subjectT2, '.');
T2name = subjectT2(1:dots(1)-1);

dots = strfind(subjectDTI, '.');
DTIname = subjectDTI(1:dots(1)-1);

clear dots

%% PID T2 to T1

switch cfg.fsl.T2interp
    case 'spline'
        %% Get transfer matrix
        
%         aa=horzcat(cfg.fsl.version,'flirt -in ',subjectT2,' -ref ',subjectT1,' -out ',T2name,'_T1Space_trilinear.nii -omat',...
%               ' T2_T1Space.mat ','-bins 256 -cost ',cfg.fsl.T2cost,' -searchrx -180 180',...
%               ' -searchry -180 180 -searchrz -180 180 -dof ',cfg.fsl.T2dof,' -interp trilinear');
%         
        cmd = sprintf('%s/flirt -in %s -ref %s -out %s_T1Space_trilinear.nii -omat T2_T1Space.mat -bins 256 -cost %s -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof %s -interp trilinear', ...
            cfg.fsl.path, ...
            subjectT2, ...
            subjectT1, ...
            T2name, ...
            cfg.fsl.T2cost, ...
            cfg.fsl.T2dof ...
            );
        
        %unix(aa,'-echo')
        unix(cmd,'-echo')
        
        %% Use nonlinear registration to use spline interpolation
        
%                  aa=horzcat(cfg.fsl.version,'applywarp --in=',subjectT2,' --ref=',subjectT1,' --out=',T2name,'_T1Space_spline.nii '...
%               ,'--premat=T2_T1Space.mat --interp=spline');
%         
        %
        cmd=sprintf('%s/applywarp --in=%s --ref=%s --out=%s_T1Space_spline.nii --premat=T2_T1Space.mat --interp=spline', ...
            cfg.fsl.path, ... 
            subjectT2, ...    
            subjectT1, ...
            T2name ...
            );
        
        % unix(aa,'-echo')
        unix(cmd,'-echo')
        
        subjectT2aligned = horzcat(T2name,'_T1Space_trilinear.nii.gz');
        T2_aligned = ft_read_mri(subjectT2aligned);
        
    case {'trilinear','nearestneighbour','sinc'}
        
        %aa=horzcat(cfg.fsl.path,'flirt -in ',subjectT2,' -ref ',subjectT1,' -out ',T2name,'_T1Space',cfg.fsl.T2interp,'.nii -omat T2_T1Space.mat ','-bins 256 -cost ',cfg.fsl.T2cost,' -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof ',cfg.fsl.T2dof,' -interp ',cfg.fsl.T2interp);
        
        cmd=sprintf('%s/flirt -in %s -ref %s -out %s_T1Space%s.nii -omat T2_T1Space.mat -bins 256 -cost %s -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof %s -interp %s', ... 
            cfg.fsl.path, ...
                subjectT2, ...
                subjectT1, ...
                T2name, ...
                cfg.fsl.T2interp, ...
                cfg.fsl.T2cost, ...
                cfg.fsl.T2dof, ...
                cfg.fsl.T2interp ...
                );
        
        %unix(aa,'-echo')
        
        subjectT2aligned=horzcat(T2name,'_T1Space',cfg.fsl.T2interp,'.nii.gz');
        T2_aligned = ft_read_mri(subjectT2aligned);
        
    otherwise
        error('unrecognized interpolation option')
        
end

%% for DTI

% aa=horzcat(cfg.fsl.path,'flirt -in ',subjectDTI,' -ref ',subjectT2aligned,' -out ',DTIname,'_T1Space_trilinear.nii -omat DTI_T1Space.mat ','-bins 256 -cost ',cfg.fsl.DTIcost,' -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof ',cfg.fsl.DTIdof,' -interp trilinear');
% unix(aa,'-echo')
% 
% delete(horzcat(DTIname,'_T1Space_trilinear.nii*'))
% 
% if length(size(DTI.anatomy))>3
%     aa=horzcat(cfg.fsl.path,'-fslsplit ',DTIname,' splitted_',DTIname,'_ -t');
%     unix(aa,'-echo')
% else
% end;
% 
% bb='';
% switch cfg.fsl.T2interp
%     case 'spline'
%         for ii=1:10
%             
%             bb=horzcat(bb,' DTI_T1_space_spline_000',num2str(ii-1));
%             aa=horzcat(cfg.fsl.path,'applywarp --in=splitted_',DTIname,'_000',num2str(ii-1),' --ref=',subjectT2aligned,' --out=DTI_T1_space_spline_000',num2str(ii-1),' --premat=DTI_T1Space.mat --interp=spline');
%             unix(aa,'-echo')
%         end;
%         for ii=11:size(DTI.anatomy,4);
%             bb=horzcat(bb,' DTI_T1_space_spline_00',num2str(ii-1));
%             aa=horzcat(cfg.fsl.path,'applywarp --in=splitted_',DTIname,'_00',num2str(ii-1),' --ref=',subjectT2aligned,' --out=DTI_T1_space_spline_00',num2str(ii-1),' --premat=DTI_T1Space.mat --interp=spline');
%             unix(aa,'-echo')
%         end;
%         
%         delete(horzcat('splitted_',DTIname,'*'))
%         
%         aa=horzcat(cfg.fsl.path,'fslmerge -t ',DTIname,'_T1_space',bb);
%         unix(aa,'-echo')
%         
%         delete(horzcat('splitted_',DTIname,'*'));
%         delete('DTI_T1_space_spline*');
%         delete(DTIname,'_T1Space_trilinear.nii*');
%         
%     case {'trilinear','nearestneighbour','sinc'}
%         for ii=1:10
%             
%             bb=horzcat(bb,' DTI_T1_space_',cfg.fsl.DTIinterp,'_000',num2str(ii-1));
%             aa=horzcat(cfg.fsl.path,'flirt -in splitted_',DTIname,'_000',num2str(ii-1),' -ref ',subjectT2aligned,' -out DTI_T1_space_',cfg.fsl.DTIinterp,'_000',num2str(ii-1),' -init DTI_T1Space.mat -applyxfm');
%             unix(aa,'-echo')
%         end;
%         for ii=11:size(DTI.anatomy,4);
%             
%             bb=horzcat(bb,' DTI_T1_space_',cfg.fsl.DTIinterp,'_00',num2str(ii-1));
%             aa=horzcat(cfg.fsl.path,'flirt -in splitted_',DTIname,'_00',num2str(ii-1),' -ref ',subjectT2aligned,' -out DTI_T1_space_',cfg.fsl.DTIinterp,'_00',num2str(ii-1),' -init DTI_T1Space.mat -applyxfm');
%             unix(aa,'-echo')
%         end;
%         
%         delete(horzcat('splitted_',DTIname,'*'))
%         aa=horzcat(cfg.fsl.path,'fslmerge -t ',DTIname,'_T1_space',bb);
%         unix(aa,'-echo')
%         
%         delete(horzcat('splitted_',DTIname,'*'));
%         delete(horzcat('DTI_T1_space_',cfg.fsl.DTIinterp,'*'));
%         delete(DTIname,'_T1Space_trilinear.nii*');
%         
%     otherwise
%         error('unrecognized interpolation option')
%         
% end
% 
% 
% %% To get the rotational part of the transformation matrix
% [cc dd]=unix(horzcat(cfg.fsl.path,'avscale DTI_T1Space.mat'),'-echo');
% 
% dummy = strfind(dd, ':');
% dummy2= strfind(dd, 'Scales');
% ee=dd(dummy(1)+1:dummy2(1)-1);
% ee=str2num(ee);
% rot_mat=ee(1:3,1:3);
% clear dummy dummy2 cc dd ee
% 
% %% rotate the diffusion directions with the rotational part of the transformation
% or_bvec=importdata(horzcat(DTIname,'.bvec'));
% rot_bvec=zeros(3,size(DTI.anatomy,4));
% for ii=1:size(DTI.anatomy,4)
%     dummy_bvec=or_bvec(:,ii);
%     rot_bvec(:,ii)=rot_mat*dummy_bvec;
% end
% 
% %write new bvec file
% fid = fopen(horzcat(DTIname,'_T1Space.bvec'), 'wt');
% fprintf(fid, '%6.8f ', rot_bvec(1,:));
% fprintf(fid,'\n');
% fprintf(fid, '%6.8f ', rot_bvec(2,:));
% fprintf(fid,'\n');
% fprintf(fid, '%6.8f ', rot_bvec(3,:));
% fclose(fid);
% 
% subjectDTIaligned=horzcat(DTIname,'_T1_space.nii.gz');
% DTI_aligned = ft_read_mri(subjectDTIaligned);
% 
% 
%%
% they are now again in memory, but realigned to each other. This is where
% a normal fieldtrip analysis pipeline would continue.
%figure; ft_plot_ortho(T1.anatomy); title('T1 image')
%figure; ft_plot_ortho(T2_aligned.anatomy); title('T2 after aligning')
% figure; ft_plot_ortho(DTI_aligned.anatomy(:,:,:,1));title('DTI after aligning')

assert(isequal(T1.transform,T2_aligned.transform),'Transformation matrix changed due to alignement.');
assert(~isequal(T2.transform,T2_aligned.anatomy),'Anatomy has not changed despite of alignement.');

T2_4p1_sinc = T2_aligned;


% ??? Error using ==> ft_read_mri at 61
% file
% '~/meeting_simbio/test_bug1826_data/test_bug1826_data/T2_T1Spacenearestneighbour.nii.gz'
% does not exist
