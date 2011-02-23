function [res] = mne_ex_evoked_grad_amp(inname,bmin,bmax,outname)
%
%   function [res] = mne_ex_evoked_grad_amp(inname,bmin,bmax,outname)
%
%   Compute the magnitude of the tangential gradient at each
%   sensor location using the planar gradiometer data and
%   optionally output the result to a fif file.
%
%   inname      The input file name. All average data sets are
%               read and processed
%   bmin,bmax   Baseline limits in seconds
%   outname     Optional output file name
%
%
%   Function returns the data which was or would have been written
%   to the file
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%

me='MNE:mne_ex_evoked_grad_amp';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin == 1
    do_baseline = false;
elseif (nargin == 3 || nargin == 4)
    do_baseline = true;
else
    error(me,'Wrong number of arguments');
end
%
%   Read the data
%
data = fiff_read_evoked_all(inname);
%
%   Figure out the planar gradiometer pairs
%
pairs = zeros(data.info.nchan,2);
npair = 0;
k = 1;
while k < data.info.nchan
    %
    %   First check the coil types
    %
    coil1 = data.info.chs(k).coil_type;
    coil2 = data.info.chs(k+1).coil_type;
    if (coil1 == coil2 && ...
            (coil1 == 2 || coil1 == 3012 || coil1 == 3013))
        one = data.info.ch_names{k};
        two = data.info.ch_names{k+1};
        lastone = one(length(one));
        lasttwo = two(length(two));
        %
        %   Then the channel names
        %
        if (strcmp(one(1:3),'MEG') && strcmp(one(1:3),'MEG'))
            if (strcmp(one(1:7),two(1:7)) && ...
                    ((lastone == '2' && lasttwo == '3') || ...
                    (lastone == '3' && lasttwo == '2')))
                npair = npair + 1;
                pairs(npair,1) = k;
                pairs(npair,2) = k+1;
                k = k + 1;
            end
        end
    end
    k = k + 1;
end

if npair == 0
    error(me,'No planar gradiometers in these data');
end
%
%   Compute the amplitudes
%
fprintf(1,'Computing the amplitudes');
if do_baseline
    fprintf(1,' (Baseline = %7.1f ... %7.1f ms)',1000*bmin,1000*bmax);
end
fprintf(1,'...');
for k = 1:length(data.evoked)
    epochs = data.evoked(k).epochs;
    %
    %  Setup baseline limits
    %
    if do_baseline
        b1 = double(data.info.sfreq*bmin - data.evoked(k).first);
        b2 = double(data.info.sfreq*bmax - data.evoked(k).first);
        if b1 < 1
            b1 = 1;
        end
        if b2 > size(epochs,2)
            b2 = size(epochs,2)
        end
    else
        b1 = 1;
        b2 = 1;
    end
    %
    %   Go through all pairs
    %
    for p = 1:npair
        one = pairs(p,1);
        two = pairs(p,2);
        if b2 > b1
            base1 = sum(epochs(one,b1:b2))/(b2-b1);
            base2 = sum(epochs(two,b1:b2))/(b2-b1);
            epochs(one,:) = sqrt((epochs(one,:)-base1).*(epochs(one, ...
                :)-base1)+(epochs(two,:)-base2).*(epochs(two,:)-base2));
        else
            epochs(one,:) = sqrt(epochs(one,:).*epochs(one, ...
                :)+epochs(two,:).*epochs(two,:));
        end
    end
    data.evoked(k).epochs = epochs;
    fprintf(1,'.');
end
fprintf(1,'[done]\n');
%
%   Compose the selection name list
%
pairs = pairs(1:npair,:);
for k = 1:npair
    ch_sel_names{k} = data.info.ch_names{pairs(k,1)};
end
%
%   Omit MEG channels but include others
%
for p = 1:data.info.nchan
    if (data.info.chs(p).kind ~= FIFF.FIFFV_MEG_CH)
        k = k + 1;
        ch_sel_names{k} = data.info.ch_names{p};
    end
end
%
%   Modify the bad channel list
%
if ~isempty(data.info.bads)
    nbad = length(data.info.bads);
    for k = 1:npair
        one = data.info.ch_names{pairs(k,1)};
        two = data.info.ch_names{pairs(k,2)};
        %
        %   If one channel of the planar gradiometer is marked bad,
        %   add the other to the bad channel list
        %
        if (~isempty(strmatch(one,data.info.bads)) && ...
                isempty(strmatch(two,data.info.bads)))
            nbad = nbad + 1;
            data.info.bads{nbad} = two;
        elseif (isempty(strmatch(one,data.info.bads)) && ...
                ~isempty(strmatch(two,data.info.bads)))
            nbad = nbad + 1;
            data.info.bads{nbad} = one;
        end
    end
end
%
%   Do the picking
%
res = fiff_pick_channels_evoked(data,ch_sel_names);
%
%   Optionally write an output file
%
if nargin == 4
    fiff_write_evoked(outname,res);
    fprintf(1,'Wrote %s\n',outname);
end
