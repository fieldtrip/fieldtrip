function [data,times] = fiff_read_raw_segment(raw,from,to,sel)
%
% [data,times] = fiff_read_raw_segment(raw,from,to,sel)
%
% Read a specific raw data segment
%
% raw    - structure returned by fiff_setup_read_raw
% from   - first sample to include. If omitted, defaults to the
%          first sample in data
% to     - last sample to include. If omitted, defaults to the last
%          sample in data
% sel    - optional channel selection vector
%
% data   - returns the data matrix (channels x samples)
% times  - returns the time values corresponding to the samples (optional)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.9  2007/11/22 05:04:25  msh
%   Fixed help text and some error checking
%
%   Revision 1.8  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.7  2006/04/21 22:11:34  msh
%   Times calculation needed casting to double
%
%   Revision 1.6  2006/04/21 22:03:43  msh
%   Fixed time limit output.
%
%   Revision 1.5  2006/04/21 16:17:48  msh
%   Added handling of CTF compensation
%
%   Revision 1.4  2006/04/21 15:11:56  msh
%   Made the debug output optional.
%
%   Revision 1.3  2006/04/21 14:43:59  msh
%   Report projection operators.
%   Some more checks in raw data reading.
%
%   Revision 1.2  2006/04/21 14:23:16  msh
%   Further improvements in raw data reading
%
%

me='MNE:fiff_read_raw_segment';

if nargin == 3
    sel = [];
elseif nargin == 2
    to  = raw.last_samp;
    sel = [];
elseif nargin == 1
    from = raw.first_samp;
    to   = raw.last_samp;
    sel  = [];
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end
%
%  Initial checks
%
from = double(from);
to   = double(to);
if from < raw.first_samp
    from = raw.first_samp;
end
if to > raw.last_samp
    to = raw.last_samp;
end
%
if from > to
    error(me,'No data in this range');
end
fprintf(1,'Reading %d ... %d  =  %9.3f ... %9.3f secs...', ...
    from,to,from/raw.info.sfreq,to/raw.info.sfreq);
%
%  Initialize the data and calibration vector
%
nchan = raw.info.nchan;
dest  = 1;
cal   = diag(raw.cals);
%
if isempty(sel)
    data = zeros(nchan,to-from+1);
    if isempty(raw.proj) && isempty(raw.comp)
        mult = [];
    else
        if isempty(raw.proj)
            mult = raw.comp*cal;
        elseif isempty(raw.comp)
            mult = raw.proj*cal;
        else
            mult = raw.proj*raw.comp*cal;
        end
    end
else
    data = zeros(length(sel),to-from+1);
    if isempty(raw.proj) && isempty(raw.comp)
        mult = [];
        cal  = diag(raw.cals(sel));
    else
        if isempty(raw.proj)
            mult = raw.comp(sel,:)*cal;
        elseif isempty(raw.comp)
            mult = raw.proj(sel,:)*cal;
        else
            mult = raw.proj(sel,:)*raw.comp*cal;
        end
    end
end
do_debug=false;
if ~isempty(cal)
    cal = sparse(cal);
end
if ~isempty(mult)
    mult = sparse(mult);
end

if isempty(fopen(raw.fid))
   fid = fopen(raw.info.filename,'rb','ieee-be');
   if (fid < 0)
      error(me,'Cannot open file %s',raw.info.filename);
   end
else
   fid = raw.fid;
end

for k = 1:length(raw.rawdir)
    this = raw.rawdir(k);
    %
    %  Do we need this buffer
    %
    if this.last > from
        if isempty(this.ent)
            %
            %  Take the easy route: skip is translated to zeros
            %
            if do_debug
                fprintf(1,'S');
            end
            doing_whole = false;
            if isempty(sel)
                one = zeros(nchan,this.nsamp);
            else
                one = zeros(length(sel),this.nsamp);
            end
        else
            tag = fiff_read_tag(fid,this.ent.pos);
            %
            %   Depending on the state of the projection and selection
            %   we proceed a little bit differently
            %
            if isempty(mult)
                if isempty(sel)
                    one = cal*double(reshape(tag.data,nchan,this.nsamp));
                else
                    one = double(reshape(tag.data,nchan,this.nsamp));
                    one = cal*one(sel,:);
                end
            else
                one = mult*double(reshape(tag.data,nchan,this.nsamp));
            end
        end
        %
        %  The picking logic is a bit complicated
        %
        if to >= this.last && from <= this.first
            %
            %    We need the whole buffer
            %
            first_pick = 1;
            last_pick  = this.nsamp;
            if do_debug
                fprintf(1,'W');
            end
        elseif from > this.first
            first_pick = from - this.first + 1;
            if to < this.last
                %
                %   Something from the middle
                %
                last_pick = this.nsamp + to - this.last;
                if do_debug
                    fprintf(1,'M');
                end
            else
                %
                %   From the middle to the end
                %
                last_pick = this.nsamp;
                if do_debug
                    fprintf(1,'E');
                end
            end
        else
            %
            %    From the beginning to the middle
            %
            first_pick = 1;
            last_pick  = to - this.first + 1;
            if do_debug
                fprintf(1,'B');
            end
        end
        %
        %   Now we are ready to pick
        %
        picksamp = last_pick - first_pick + 1;
        if picksamp > 0
            data(:,dest:dest+picksamp-1) = one(:,first_pick:last_pick);
            dest = dest + picksamp;
        end
    end
    %
    %   Done?
    %
    if this.last >= to
        fprintf(1,' [done]\n');
        break;
    end
end

fclose(fid);

if nargout == 2
    times = [ from:to ];
    times = double(times)/raw.info.sfreq;
end

return;

end
