function nmt_image(op,varargin)
% analogous to spm's spm_image -- used to update coords when clicking
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
            fg  = spm_figure('Findwin','Graphics');
            
            %%
            posmrimm = spm_orthviews('pos')';
 %           set(st.nmt.gui.megp,'String',sprintf('%.1f %.1f %.1f',posmrimm));
           
            if(prod(size(st.vols{1}.blobs{1}.vol)) == length(st.nmt.pos))
                % this is the case for official Fieldtrip output based on cube/rectangular VOIs 
                blobidx = nmt_transform_coord(inv(st.vols{1}.blobs{1}.mat),posmrimm);
                blobidx = round(blobidx);
                blobdim = size(st.vols{1}.blobs{1}.vol);
                if(all(blobidx <= blobdim & blobidx > 0))
                    %                blobintensity = st.vols{1}.blobs{1}.vol(blobidx(1),blobidx(2),blobidx(3));
                    st.nmt.cfg.vox_idx  = sub2ind(size(st.vols{1}.blobs{1}.vol),blobidx(1),blobidx(2),blobidx(3));
                    blobintensity = st.nmt.fun{1}(st.nmt.cfg.vox_idx);
                else % if out of blob's bounds
                    blobintensity = NaN;
                    st.nmt.cfg.vox_idx = NaN;
                end
            else
                % otherwise, this is an arbitrary coordinate list (e.g., Nutmeg-style)
                [st.nmt.cfg.vox_idx, d] = dsearchn(st.nmt.pos, posmrimm);
                if (d < st.vols{1}.blobs{1}.mat(1)/sqrt(2))
                    % st.vols{1}.blobs{1}.mat is the functional voxelsize
                    % if selected voxel is farther than sqrt(1/2)*voxelsize from
                    % nearest functional voxel, then we must be outside the VOI
                    blobintensity = st.nmt.fun{1}(st.nmt.cfg.vox_idx);
                else
                    blobintensity = NaN;
                end
            end
            
            set(st.nmt.gui.beamin,'String',sprintf('%+g',blobintensity));
            
            if(~isempty(st.nmt.cfg.atlas))
                atlas_labels = atlas_lookup(st.nmt.cfg.atlas,posmrimm,'inputcoord','mni','queryrange',3);
                set(st.nmt.gui.mnilabel,'String',atlas_labels);
            end

            
            if(size(st.nmt.fun{1},2)>1)
                nmt_timeselect;
                
                st.nmt.cfg.voxinside_idx=find(st.nmt.cfg.vox_idx==st.nmt.cfg.inside_idx);
                if(isempty(st.nmt.cfg.voxinside_idx)), st.nmt.cfg.voxinside_idx=NaN,end
%flash_statori_plotori(st.nmt.cfg.voxinside_idx);
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


