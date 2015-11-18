function nmt_image(op,varargin)
% analagous to spm's spm_image -- used to update coords when clicking
% around MRI viewer

global st

% reset MRI axis positions to proper location, in case they moved (e.g.,
% due to manipulation by SPM controls)
for ii=1:3
    set(st.vols{1}.ax{ii}.ax,'Position',[st.nmt.gui.mriaxpos(:,ii)]);
end


feval(st.nmt.spm_refresh_handle);

if(~exist('op','var') || ~ischar(op))
    op = 'shopos';
end

switch(op)
    case 'shopos'
        if isfield(st,'mp'),
            fg  = spm_figure('Findwin','Graphics');
%             posmrimm = spm_orthviews('pos')';
%             set(st.nmt.gui.megp,'String',sprintf('%.1f %.1f %.1f',posmrimm));
           
             
%             blobidx = nmt_transform_coord(inv(st.vols{1}.blobs{1}.mat),posmrimm);
%             blobidx = round(blobidx);
%             blobdim = size(st.vols{1}.blobs{1}.vol);
%             if(all(blobidx <= blobdim & blobidx > 0))
%                 blobintensity = st.vols{1}.blobs{1}.vol(blobidx(1),blobidx(2),blobidx(3));
%                 st.nmt.cfg.vox_idx  = sub2ind(size(st.vols{1}.blobs{1}.vol),blobidx(1),blobidx(2),blobidx(3));
%             else % if out of blob's bounds
%                 blobintensity = NaN;
%                 st.nmt.cfg.vox_idx = NaN;
%             end
%             set(st.nmt.gui.beamin,'String',sprintf('%+g',blobintensity));
            
            
            if(size(st.nmt.fun,2)>1)
                nmt_timeselect;
            end
            
%             [~,mriname]=fileparts(st.vols{1}.fname);
%             if isfield(coreg,'norm_mripath')
%                 dbout=dbstack;
%                 if length(dbout)>1 & ( strcmp(dbout(2).name,'paint_activation') )
%                     posmni = nut_mri2mni(posmrimm,coreg,0);
%                 else
%                     posmni = nut_mri2mni(posmrimm,coreg,1);
%                 end
%                 set(st.mnip,'String',sprintf('%.1f %.1f %.1f',posmni));
%             elseif strcmp(mriname(1),'w')
%                 posmni = posmrimm;
%                 set(st.mnip,'String',sprintf('%.1f %.1f %.1f',posmni));
%             end
%             
%             if isfield(coreg,'norm_mripath') || strcmp(mriname(1),'w')
%                 if(isfield(rivets,'TalDB') && exist('posmni','var'))
%                     ind = rivets.TalDB.data(dsearchn(rivets.TalDB.coords,nut_mni2tal(posmni)),:);
%                     if(ind==0)
%                         set(st.mnilabel,'String','');
%                     else
%                         set(st.mnilabel,'String',[rivets.TalDB.labels{:,ind}]);
%                     end
%                 else
%                     set(st.mnip,'String','');
%                 end
%             end
%                 
%                 % set sliders
%                 if ndefaults.sliders
%                 posmrivx = nut_mm2voxels(posmrimm);
%                 
%                 mat = st.vols{1}.premul*st.vols{1}.mat;
%                 R=spm_imatrix(st.vols{1}.mat);
%                 R = spm_matrix([0 0 0 R(4:6)]);
%                 R = R(1:3,1:3);
%                 dim = (st.vols{1}.dim*R)';
%                 
%                 for i=1:3
%                     set(st.slider{i},'Value',posmrivx(i)/dim(i));
%                 end
%             end
% 
%             if(~isempty(beam))
%                 feval(rivets.ts_refresh_handle);
%             end
%             
%             if isfield(beam,'corr') && isfield(beam.corr,'behavdata')      % Adrian added this for the future new FCM toolbox
%                 try, fcm_showcorr; end
%             end
            
        else
            st.Callback = ';';
        end
        return
    case 'setposmeg'
        warning('TODO: not yet fully implemented');
        if isfield(st.nmt.gui,'megp'),
            fg = spm_figure('Findwin','Graphics');
            if(any(findobj(fg) == st.nmt.gui.megp))
                pos = sscanf(get(st.nmt.gui.megp,'String'), '%g %g %g',[1 3]);
                if(length(pos)==3)
                    pos = nut_meg2mri(pos);
                else
                    pos = spm_orthviews('pos');
                end
                nmt_repos(pos);
            end
        end
        return
    case 'setposmni'
        warning('TODO: not yet fully implemented');
        if isfield(st.nmt.gui,'mnip'),
            fg = spm_figure('Findwin','Graphics');
            if(any(findobj(fg) == st.nmt.mnip))
                pos = sscanf(get(st.nmt.mnip,'String'), '%g %g %g',[1 3]);
                if(length(pos)==3)
                    pos = nut_mni2mri(pos,coreg);
                else
                    pos = spm_orthviews('pos');
                end
                nmt_repos(pos);
            end
        end
        return
    otherwise
        error('you''re doing something wrong.');
end


