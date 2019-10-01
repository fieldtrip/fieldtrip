% movav() - Perform a moving average of data indexed by xvals.
%           Supports use of a moving non-rectangular window.
%           Can be used to resample a data matrix to any size 
%           (see xadv NOTE below) and to regularize sampling of 
%           irregularly sampled data.
% Usage:
%  >> [outdata,outx] = movav(data,xvals,xwidth,xadv,firstx,lastx,xwin,nonorm);
%
% Input:
%   data   = input data (chans,frames)
%
% Optional inputs:
%   xvals  = increasing x-values for data frames (columnsa). The default 
%            [1:frames] is fastest {def|[]|0 -> 1:frames}
%   xwidth = smoothing-window width in xvals units {def|0->(hix-lox)/4}
%   xadv   = window step size in xvals units. NOTE: To reduce yyy frames 
%            to about xxx, xadv needs to be near yyy/xxx {default|0 -> 1}
%   firstx = low xval of first averaging window {def|[] -> min xvals}
%   lastx  = high xval of last averaging window {def|[] -> max xvals}
%   xwin   = vector of window values {def|0 -> ones() = square window}
%            May be long. NOTE: linear interp. is NOT used between values.
%            Example: gauss(1001,2) ->  [0.018 ... 1.0 ... 0.018]
%   nonorm = [1|0] If non-zero, do not normalize the moving sum. If
%            all y values are 1s. this creates a moving histogram. 
%            Ex: >> [oy,ox] = movav(ones(size(x)),x,xwd,xadv,[],[],0,1);
%            returns a moving histogram of xvals {default: 0}
% Outputs:
%   outdata = smoothed output data (chans,outframes)
%   outx    = xval midpoints of successive output windows
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 10-25-97 

% Copyright (C) 10-25-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 3-20-98 fixed bug in multi-channel windowed averaging -sm
% 6-10-98 changed mean() and sum() to nanmean() and nansum() -sm
% 2-16-99 tested for stat toolbox functions nanmean() and nansum() -sm
% 9-03-01 fixed gauss() example -sm
% 01-25-02 reformated help & licenses -ad 

function [outdata,outx] = movav(data,xvals,xwidth,xadv,firstx,lastx,xwin,nonorm)

MAXPRINT = 1;     % max outframe numbers to print on tty
NEARZERO = 1e-22; %
DEFAULT_XADV = 1; % default xvals window step advance
verbose = 0;      % If 1, output process info

nanexist = 0;  
if nargin<1
   help movav
   return
else
  [chans,frames]=size(data);
end
if chans>1 & frames == 1,
  data   = data';   % make row vector
  tmp    = chans;
  chans  = frames;
  frames = tmp;
end
if frames < 4
  error('data are too short');
  return
end

flag_fastave = 0;
if nargin<2 | isempty(xvals) | (numel(xvals)==1 & xvals == 0)
  xvals = 1:frames; % flag default xvals
  flag_fastave = 0; % TURNED OFF THIS FEATURE - LEADS TO ?? BUG AT ABOUT 287
end                 % -sm 3/6/07
if size(xvals,1)>1 & size(xvals,2)>1
  error('xvals must be a vector');
end
xvals = xvals(:)'; % make xvals a row vector

if frames ~= length(xvals) 
    error('lengths of xvals and data not equal');
end

if nargin < 8 | isempty(nonorm)
  nonorm = 0;  % default -> return moving mean
end
if abs(nonorm) > NEARZERO
   nonorm = 1;
end

if nargin < 7 | isempty(xwin)
  xwin = 0;
end

if nargin < 6 | isempty(lastx)
  lastx = [];
end
if isempty(lastx),
  if flag_fastave
    lastx = frames;
  else
    lastx = max(xvals);
  end
end

if nargin<5 | isempty(firstx)
  firstx = [];
end
if isempty(firstx),
  if flag_fastave
    firstx = 1;
  else
    firstx = min(xvals);
  end
end

if nargin<4 | isempty(xadv)
  xadv = 0;
end
if isempty(xadv) | xadv == 0,
  xadv = DEFAULT_XADV;
end

if nargin<3 | isempty(xwidth) | xwidth==0
  xwidth = (lastx-firstx)/4;  % DEFAULT XWIDTH
end

wlen = 1;  % default;
if flag_fastave==0
  if length(xwin)==1 & (xwin~=0) & (xwin~=1),  % should be a vector or 0
    error('xwin not vector or 0');
  elseif size(xwin,1)>1 & size(xwin,2)>1 % not a matrix
    error('xwin cannot be a matrix'); 
  end
  if size(xwin,1)>1
    xwin = xwin';   % make row vector
  end

  if xwin~=0
    wlen = length(xwin);
  end
end

%outframes = floor(0.99999+((lastx-firstx+xadv)-xwidth)/xadv);
outframes = floor(((lastx-firstx+xadv+1)-xwidth)/xadv);
if verbose
  fprintf('movav() will output %d frames.\n',outframes);
end
if outframes < 1,
   outframes = 1;
end
outdata = zeros(chans,outframes);
outx = zeros(1,outframes);
outxval = firstx+xwidth/2;
%
%%%%%%%%%%%%%%%%%%%%%% Print header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
  fprintf('Performing moving averaging:\n')
  fprintf('Output will be %d chans by %d frames',chans,outframes);
  if wlen>1,
    fprintf(' using the specified width-%d window\n',wlen);
  else
    fprintf(' using a width-%d square window\n',xwidth);
  end
  fprintf(' and a window advance of %g\n',xadv);
end
%
%%%%%%%%%%%%%%%%%%% Perform averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
lox = firstx;
i = 0; % flag_fastave default
for f=1:outframes
   hix = lox+xwidth;
   outx(1,f)=outxval;
   outxval = outxval + xadv;
   if flag_fastave == 0 
      i = find(xvals>=lox & xvals < hix);
   end
   if length(i)==0,
      if f>1,
       outdata(:,f) = outdata(:,f-1); % If no data, replicate
      else
       outdata(:,f) = zeros(chans,1); %  or else output zeros
      end
   elseif length(xwin)==1,
      if flag_fastave > 0
          outdata(:,f) = nan_mean(data(:,round(lox):round(hix))')'; 
          nix = length([round(lox):round(hix)]);
      else
          outdata(:,f) = nan_mean(data(:,i)')'; % Else average
          nix = length(i);
      end
      if nonorm & nix % undo division by number of elements summed
          outdata(:,f) = outdata(:,f)*nix;
      end
%
%%%%%%%%%%%%%%%%% Windowed averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   else % length(xwin) > 1
       wadv=(hix-lox)/wlen;
       ix = floor((xvals(i)-lox)/wadv)+1; % AG fix 3/6/07
       if length(xwin)>1
          sumx = sum(xwin(ix));
       else
          sumx=1;
       end

       % AG fix 3/6/7
       outdata(:,f) = nan_sum((((ones(chans,1)*xwin(ix)).*data(:,i)))')';
       if abs(sumx) > NEARZERO & nonorm == 0 
          outdata(:,f) = outdata(:,f)/sumx;
       end
   end
   lox = lox+xadv;
   if (outframes<MAXPRINT) 
      fprintf('%d ',f);
   end
end
if verbose,
  fprintf('\n');
end

%
%%%%%%%%%%%%%%%%%%%%%%% function nan_mean() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% nan_mean() - take the column means of a matrix, ignoring NaN values
%
function out = nan_mean(in)

   nans = find(isnan(in));
   in(nans) = 0;
   sums = sum(in);
   nonnans = ones(size(in));
   nonnans(nans) = 0;
   nonnans = sum(nonnans,1);
   nononnans = find(nonnans==0);
   nonnans(nononnans) = 1;
   out = sum(in,1)./nonnans;
   out(nononnans) = NaN;

%
%%%%%%%%%%%%%%%%%%%%%%% function nan_sum() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% nan_sum() - take the column sums of a matrix, ignoring NaN values
%
function out = nan_sum(in)

   nans = find(isnan(in));
   in(nans) = 0;
   out = sum(in,1);

   nonnans = ones(size(in));
   nonnans(nans) = 0;
   nonnans = sum(nonnans,1);
   nononnans = find(nonnans==0);
   out(nononnans) = NaN;
