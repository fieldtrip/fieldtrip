function nmt_tfplot(axh,twin,fwin,tf)
% nutmegtrip uses this helper function to plot time-frequency data

% reformat beam.bands and corresponding tf info to include gaps
t = unique(twin(:));
f = fwin';
f = f(:);

for ii=1:(length(f)-1)
    fwinfin(ii,:) = [f(ii) f(ii+1)];
end
ffin = fwinfin';
ffin = ffin(:);

[~,fwinfin2fwin]=intersect(fwinfin,fwin,'rows');

z = nan(size(ffin,1),size(tf,1));

for ii=1:size(fwin,2)
    z(2*fwinfin2fwin(ii)-1,:) = tf(:,ii);
    z(2*fwinfin2fwin(ii),:) = tf(:,ii);
end

global st
st.nmt.gui.h_tf = mesh(axh,t,ffin,z);
view(axh,2); % '2-D view' of spectrogram
%set(st.nmt.gui.h_tf,'LineStyle','none'); % useful if plot made with 'surf'
axis(axh,'tight');
set(axh,'YScale','log');

if(verLessThan('matlab','8.4')) % necessary to preserve colormap on functional image for Matlab R2014a and earlier
    for ii=1:3
        freezeColors(st.vols{1}.ax{ii}.ax);
    end
end

caxis(axh,[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max*33/32]);
colormap(axh,[jet(128); 1 1 1]);

% xlabel(axh,'Time');
ylabel(axh,'Frequency (Hz)');
