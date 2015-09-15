function nmt_image(op,varargin)
% analagous to spm's spm_image -- used to update coords when clicking
% around MRI viewer

global st

% reset MRI axis positions to proper location, in case they moved (e.g.,
% due to manipulation by SPM controls)
for ii=1:3
    set(st.vols{1}.ax{ii}.ax,'Position',[st.nmt.mriaxpos(:,ii)]);
end


feval(st.nmt.spm_refresh_handle);

if(~exist('op','var') || ~ischar(op))
    op = 'shopos';
end

switch(op)
    case 'shopos'
        if isfield(st,'mp'),
            fg  = spm_figure('Findwin','Graphics');
            posmrimm = spm_orthviews('pos')';
            set(st.nmt.megp,'String',sprintf('%.1f %.1f %.1f',posmrimm));
            
            blobidx = nmt_transform_coord(inv(st.vols{1}.blobs{1}.mat),posmrimm);
            blobidx = round(blobidx);
            blobdim = size(st.vols{1}.blobs{1}.vol);
            if(all(blobidx <= blobdim & blobidx > 0))
                blobintensity = st.vols{1}.blobs{1}.vol(blobidx(1),blobidx(2),blobidx(3));
                st.nmt.voxind  = sub2ind(size(st.vols{1}.blobs{1}.vol),blobidx(1),blobidx(2),blobidx(3));
            else % if out of blob's bounds
                blobintensity = NaN;
                st.nmt.voxind = NaN;
            end
            set(st.nmt.beamin,'String',sprintf('%g',blobintensity));
            
            
            if(size(st.nmt.fun,2)>1)
                nmt_reposition;
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
        if isfield(st,'megp'),
            fg = spm_figure('Findwin','Graphics');
            if(any(findobj(fg) == st.nmt.megp))
                pos = sscanf(get(st.nmt.megp,'String'), '%g %g %g',[1 3]);
                if(length(pos)==3)
                    pos = nut_meg2mri(pos);
                else
                    pos = spm_orthviews('pos');
                end
                nmt_reposition(pos);
            end
        end
        return
    case 'setposmni'
        warning('TODO: not yet fully implemented');
        if isfield(st,'mnip'),
            fg = spm_figure('Findwin','Graphics');
            if(any(findobj(fg) == st.nmt.mnip))
                pos = sscanf(get(st.nmt.mnip,'String'), '%g %g %g',[1 3]);
                if(length(pos)==3)
                    pos = nut_mni2mri(pos,coreg);
                else
                    pos = spm_orthviews('pos');
                end
                nmt_reposition(pos);
            end
        end
        return
    otherwise
        error('you''re doing something wrong.');
end


function nmt_reposition(op)
global st

if(~exist('op','var'))
    op = '';
end

switch(op)
    case 'ts'
        % if this was an interactive button-press in time series
        currpt = get(st.nmt.ax_ts,'CurrentPoint');
        currpt = currpt(1,1);
        
        if strcmpi(get(gcbf,'SelectionType'),'alt') % right-click
            st.nmt.cfg.time(2) = dsearchn(st.nmt.time',currpt);
        else % left-click
            st.nmt.cfg.time(1) = dsearchn(st.nmt.time',currpt);
            st.nmt.cfg.time(2) = dsearchn(st.nmt.time',currpt);
        end
    otherwise
        % otherwise, time interval was already specified, so nothing to do
end

if(st.nmt.cfg.time(2) < st.nmt.cfg.time(1))
    % enforce that time1 < time2
    st.nmt.cfg.time = st.nmt.cfg.time([2 1]);
end



if(isfinite(st.nmt.voxind))
    plot(st.nmt.ax_ts,st.nmt.time,st.nmt.fun(st.nmt.voxind,:));
    grid(st.nmt.ax_ts,'on');
else
    plot(st.nmt.ax_ts,st.nmt.time,nan(length(st.nmt.time),1));
end

% set y-axis range based on whole volume range (single time point), or
% selected timeslice range
if(st.nmt.cfg.time(1)==st.nmt.cfg.time(2))
    ymin = min(st.nmt.fun(:));
    ymax = max(st.nmt.fun(:));
else
    ymin =  min(min(st.nmt.fun(:,st.nmt.cfg.time(1):st.nmt.cfg.time(2))));
    ymax =  max(max(st.nmt.fun(:,st.nmt.cfg.time(1):st.nmt.cfg.time(2))));
end
set(st.nmt.ax_ts,'XLim',st.nmt.time([1 end]),'YLim',[ymin ymax]);
xlabel(st.nmt.ax_ts,['Time (s)']);


%% plot vertical line indicating selected time point
axisvalues = axis(st.nmt.ax_ts);
line_max=axisvalues(4);
line_min=axisvalues(3);
ts_xhair_color = 'red';
ts_xhair_width = 1;

ylim=get(st.nmt.ax_ts,'YLim');
axes(st.nmt.ax_ts);
h=patch([st.nmt.time(st.nmt.cfg.time(1)) st.nmt.time(st.nmt.cfg.time(2)) st.nmt.time(st.nmt.cfg.time(2)) st.nmt.time(st.nmt.cfg.time(1))]',[ylim(1) ylim(1) ylim(2) ylim(2)]',[1 0.4 0.4],'EdgeColor','red');

set(st.nmt.ax_ts,'ButtonDownFcn',@nmt_repos_start);
nmt_spm8_plot



function nmt_repos_start(varargin)
global st
% if(st.nmt.ax_ts == gcbo) % if ts axes were clicked (i.e., not MRI)
    set(gcbf,'windowbuttonmotionfcn',@nmt_repos_move, 'windowbuttonupfcn',@nmt_repos_end);
    nmt_reposition('ts');
% end
%_______________________________________________________________________
%_______________________________________________________________________
function nmt_repos_move(varargin)
nmt_reposition('ts');
%_______________________________________________________________________
%_______________________________________________________________________
function nmt_repos_end(varargin)
set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
