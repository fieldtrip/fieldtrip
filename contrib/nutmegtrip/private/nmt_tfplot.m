function h_tf = nmt_tfplot(axh,twin,fwin,tf,zlim,ButtonDownFcn)
% nutmegtrip uses this helper function to plot time-frequency data

% reformat beam.bands and corresponding tf info to include gaps
t = unique(twin(:));
f = fwin';
f = f(:);

for ii=1:(length(f)-1)
    fwinfin(ii,:) = [f(ii) f(ii+1)];
end
badrows = find(fwinfin(:,1)==fwinfin(:,2));
fwinfin(badrows,:) = [];
ffin = fwinfin';
ffin = ffin(:);

[~,fwinfin2fwin]=intersect(fwinfin,fwin,'rows');

z = nan(size(ffin,1),size(tf,1));

for ii=1:size(fwin,1)
    z(2*fwinfin2fwin(ii)-1,:) = tf(:,ii);
    z(2*fwinfin2fwin(ii),:) = tf(:,ii);
end

global st

h_tf = surf(axh,t,ffin,z);
view(axh,2); % '2-D view' of spectrogram
set(h_tf,'LineStyle','none'); % useful if plot made with 'surf'
set(h_tf,'ButtonDownFcn',ButtonDownFcn); % click on TF plot triggers CallbackFcn

%set(h_tf,'LineWidth',2.5); % if mesh function used to plot, helps prevent faint lines around each TF datapoint

% limit labels to defined frequencies
ytick = unique(f);
if(length(ytick)<20) % but only if there aren't too many frequency bands!
    set(axh,'YTick',unique(f));
end
axis(axh,'tight');
grid(axh,'off');
set(axh,'YScale','log');


if(verLessThan('matlab','8.4')) % necessary to preserve colormap on functional image for Matlab R2014a and earlier
    freezeColors(st.vols{1}.blobs{1}.cbar);
    for ii=1:3
        freezeColors(st.vols{1}.ax{ii}.ax);
    end
end

%caxis(axh,[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max*17/16]);
caxis(axh,[zlim(1) zlim(2)]);
colormap(axh,st.nmt.cfg.colormap);
set(axh,'Color','none'); % make blank bits transparent

% xlabel(axh,'Time');
ylabel(axh,'Frequency (Hz)');
