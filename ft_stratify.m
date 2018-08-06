function [varargout] = ft_stratify(cfg, varargin)

% FT_STRATIFY tries to reduce the variance in a specific feature in the data
% that is not related to an effect in two or multiple conditions, but where
% that feature may confound the analysis. Stratification is implemented by
% randomly removing elements from the data, making the distribution of the
% data equal on that feature.
%
% Use as
%   [output]          = ft_stratify(cfg, input1, input2, ...), or
%   [output, binaxis] = ft_stratify(cfg, input1, input2, ...)
%
% For the histogram and the split method, each input is a Nchan X Nobs
% matrix. The output is a cell-array with in each cell the same data as in
% the corresponding input, except that the observations that should be
% removed are marked with a NaN.
%
% For the equatespike method, each input is a Ntrials X 1 cell-array. Each
% trial should contain the spike firing moments (i.e. a logical Nchans X
% Nsamples matrix). The output is a cell-array with in each cell the same
% data as in the corresponding input, except that spike numbers have been
% equated in each trial and channel.
%
% The configuration should contain
%   cfg.method      = 'histogram'
%                     'splithilo'
%                     'splitlohi'
%                     'splitlolo'
%                     'splithihi'
%                     'equatespike'
%
% The following options apply only to histogram and split methods.
%   cfg.equalbinavg = 'yes'
%   cfg.numbin      = 10
%   cfg.numiter     = 2000
%
% The following options apply only to the equatespike method.
%   cfg.pairtrials  = 'spikesort', 'linkage' or 'no' (default = 'spikesort')
%   cfg.channel     = 'all' or list with indices ( default = 'all')
%
% See also FT_FREQSTATISTICS, FT_TIMELOCKSTATISTICS, FT_SOURCESTATISTICS

% Copyright (C) 2007 Jan-Mathijs Schoffelen & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% input1 and input2 are the to be stratified with respect to each other
% dimensionality of input1 (2) = chan x rpt. If nchan>1, do a "double"
% stratification

% set the defaults
cfg.method       = ft_getopt(cfg, 'method', 'histogram');
cfg.equalbinavg  = ft_getopt(cfg, 'equalbinavg', 'yes');
cfg.numbin       = ft_getopt(cfg, 'numbin', 10);
cfg.numiter      = ft_getopt(cfg, 'numiter', 2000);
cfg.binedges     = ft_getopt(cfg, 'binedges', []);
cfg.pairtrials   = ft_getopt(cfg, 'pairtrials', 'spikesort');
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.linkage      = ft_getopt(cfg, 'linkage', 'complete');      % 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'

% the input data is a cell-array containing matrices for each condition
input = varargin;

if size(input{1},1)~=size(input{2},1)
  ft_error('the number of channels should be the same');
end
if size(input{1},1)~=1 && strcmp(cfg.equalbinavg, 'yes')
  ft_error('combining equalising bin-averages and simultaneous stratification of two channels is impossible');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.method, 'histogram')
  
  nchan = size(input{1},1);
  ncond = length(input);
  % if nchan~=2, ft_error('number of trials ~= 2, do not know how to stratify'); end
  if isfield(cfg,'cmbindx')
    [output, b] = stratify2(cfg, input{1}, input{2});
    varargout{1} = output{1};
    varargout{2} = output{2};
    varargout{3} = b;
  else
    if ~isempty(cfg.binedges)
      fprintf('using the bin edges for the histogram as defined in the cfg\n');
      if iscell(cfg.binedges)
        cfg.numbin = cellfun(@numel, cfg.binedges)-1;
      else
        cfg.numbin = numel(cfg.binedges)-1;
      end
    end
    
    %------prepare some stuff
    if numel(cfg.numbin) ~= nchan
      cfg.numbin = repmat(cfg.numbin(1), [1 nchan]);
    end
    for j = 1:nchan
      tmp  = [];
      for cndlop = 1:ncond
        tmp = cat(2, tmp, input{cndlop}(j,:));
      end % concatenate input-data
      
      % create a 'common binspace'
      [ndum,x] = hist(tmp, cfg.numbin(j));
      dx    = diff(x);
      if ~isempty(cfg.binedges)
        if iscell(cfg.binedges)
          xc = cfg.binedges{j};
        else
          xc = cfg.binedges;
        end
      else
        xc    = [-inf x(1:end-1)+0.5*dx(1) inf];
      end
      
      %make histograms and compute the 'target histogram', which
      %will equalize the conditions while throwing away as few
      %repetitions as possible
      tmp = [];
      for cndlop = 1:ncond
        [n{cndlop}(j,1:numel(xc)), b{cndlop}(j,:)] = histc(input{cndlop}(j,:), xc);
        tmp = [tmp; n{cndlop}(j,:)];
      end
      binaxis(j,1:(cfg.numbin(j)+1)) = xc;
    end
    
    %------index the trials
    %------create n-d histo
    linearhisto = zeros(ncond, prod(cfg.numbin));
    for cndlop = 1:ncond
      tmpb = zeros(1, size(b{cndlop},2));
      for j = 1:nchan
        if j == 1
          tmpb = tmpb + (b{cndlop}(j,:)).*prod(cfg.numbin(1:(j-1)));
        else
          tmpb = tmpb + (b{cndlop}(j,:)-1).*prod(cfg.numbin(1:(j-1)));
        end
      end
      b{cndlop}             = tmpb;
      for binlop = 1:size(linearhisto,2)
        linearhisto(cndlop,binlop) = sum(tmpb==binlop);
      end
    end
    
    %------find intersection
    numok  = min(linearhisto,[],1);
    nummax = max(linearhisto,[],1);
    for j = 1:length(numok)
      minind{j} = find(linearhisto(:,j)==numok(j));
    end
    
    %------do stratification the easy or the hard way
    if strcmp(cfg.equalbinavg, 'yes')
      %---this is the hard way
      if nchan>1, ft_error('the option equalbinavg only works for one channel input at the moment'); end
      
      %---first allocate some memory
      for cndlop = 1:ncond
        sel{cndlop}     = zeros(1,size(b{cndlop},2));
      end
      
      for binlop = 1:length(numok)
        if numok(binlop)>0
          tmpmatmin    = nan(ncond,nummax(binlop));
          tmpmatmax    = nan(ncond,nummax(binlop));
          tmpmatminind = nan(ncond,nummax(binlop));
          tmpmatmaxind = nan(ncond,nummax(binlop));
          for cndlop = 1:ncond
            tmpsel          = find(b{cndlop}==binlop);
            tmpdat          = input{cndlop}(tmpsel);
            [tmpsrt,tmpind] = sort(tmpdat);
            tmpmatmin(   cndlop,1:linearhisto(cndlop,binlop)        ) = tmpsrt;
            tmpmatmax(   cndlop,end-linearhisto(cndlop,binlop)+1:end) = tmpsrt;
            tmpmatminind(cndlop,1:linearhisto(cndlop,binlop)        ) = tmpind;
            tmpmatmaxind(cndlop,end-linearhisto(cndlop,binlop)+1:end) = tmpind;
          end
          refavg = nanmean(tmpmatmin,2);
          refmin = mean(tmpmatmin(:,1:numok(binlop)        ),2);
          refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
          
          %---determine the amount of trials in this bin, so that for all conditions
          %it lies between the lowest possible and highest possible realisation
          ok     = 0; lowok = 0; hiok  = 0; cnt   = 0;
          offset = zeros(1,ncond);
          while ok==0
            if numok(binlop)>0
              [tmpref,refind] = min(refavg(minind{binlop}));
              if any(refmin - tmpref > 0)
                numok(binlop)                          = numok(binlop) - 1;
                offset(minind{binlop}(refind))         = offset(minind{binlop}(refind)) + 1; %correction term
                tmpmatmin(minind{binlop}(refind),:)    = [   tmpmatmin(minind{binlop}(refind),2:end) nan];
                %tmpmatminind(minind{binlop}(refind),:) = [tmpmatminind(minind{binlop}(refind),2:end) nan];
                warning off;
                refavg = nanmean(tmpmatmin,2);
                refmin = mean(tmpmatmin(:,1:numok(binlop)),2);
                refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
                warning on;
              else
                lowok = 1;
              end
              [tmpref,refind] = min(refavg(minind{binlop}));
              if any(refmax - tmpref < 0)
                numok(binlop)                          = numok(binlop) - 1;
                tmpmatmax(minind{binlop}(refind),:)    = [nan    tmpmatmax(minind{binlop}(refind),1:end-1)];
                %tmpmatmaxind(minind{binlop}(refind),:) = [nan tmpmatmaxind(minind{binlop}(refind),1:end-1)];
                warning off;
                refavg = nanmean(tmpmatmax,2);
                refmin = mean(tmpmatmin(:,1:numok(binlop)),2);
                refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
                warning on;
              else
                hiok = 1;
              end
            end
            if lowok==1 && hiok==1, ok = 1; end
            if cnt==100,            ok = 1; end
            cnt = cnt+1;
          end
          
          %---this is now the average that should be approximated
          [tmpref,refind] = min(refavg(minind{binlop}));
          
          if numok(binlop)>0
            for cndlop = 1:ncond
              pnt     = tmpmatmin(cndlop, 1:linearhisto(cndlop,binlop)) - tmpref;
              nrow    = length(pnt)-numok(binlop)+1;
              pntmat  = repmat(pnt,[nrow 1]);
              % get a good initial guess
              cpnt    = conv2([1],ones(1,numok(binlop))./numok(binlop),pnt,'same');
              indc    = nearest(cpnt,0);
              indvec  = [indc-floor(numok(binlop)/2):indc+ceil(numok(binlop)/2)-2];
              if length(indvec)<=1,       indvec = [indc indc];                    end
              if indvec(1)<1,             indvec = indvec-indvec(1)+1;             end
              if indvec(end)>length(cpnt),  indvec = indvec-indvec(end)+length(cpnt);  end
              tmpmat  = zeros(nrow,length(pnt));
              tmpmat(:, indvec                       ) = 1;
              if length(unique(indvec))>1 || size(pntmat,2)>nrow
                tmpmat(:, setdiff(1:length(pnt),indvec)) = eye(nrow);
              else
                tmpmat = eye(nrow);
              end
              %tmpmat  = [ones(nrow,numok(binlop)-1) eye(nrow)];
              %tmpmat  = tmpmat(:,randperm(size(tmpmat,2)));
              if cndlop~=minind{binlop}(refind)
                m      = nan(1,100);
                for rndlop = 1:100
                  if rndlop<=12 || sum(diff(m(rndlop-11:rndlop-1))==0)<10
                    dif = abs(sum(pntmat.*tmpmat,2)./numok(binlop));
                    [m(rndlop),ind] = min(dif);
                    tmpvec           = tmpmat(ind,:);
                    tmpmat           = repmat(tmpvec,[nrow 1]);
                    indone           = find(tmpmat(1,:));
                    tmpsel           = randperm(length(indone));
                    tmpsel           = indone(tmpsel(1));
                    tmpmat(:,[tmpsel find(tmpmat(1,:)==0)]) = eye(nrow);
                  else
                    %do nothing and go on
                    break
                  end
                end
              else
                tmpvec = [ones(1,numok(binlop)) zeros(1,nrow-1)];
              end
              tmpsel = find(b{cndlop}==binlop);
              sel{cndlop}(tmpsel(tmpmatminind(cndlop,find(tmpvec)))) = 1;
            end
          end
        end
      end
      
    else
      %------stratify the easy way
      for cndlop = 1:ncond
        sel{cndlop} = zeros(1,size(b{cndlop},2));
        tmphisto    = linearhisto(cndlop,:);
        for binlop = 1:size(tmphisto,2)
          tmpsel = find(b{cndlop}==binlop);
          tmpsel = tmpsel(randperm(length(tmpsel)));
          tmpsel = tmpsel(1:numok(binlop));
          sel{cndlop}(tmpsel) = 1;
        end
      end
    end
    
    %------create output
    output = input;
    for cndlop = 1:ncond
      output{cndlop}(:,find(sel{cndlop}==0))=nan;
    end
    varargout{1} = output;
    varargout{2} = binaxis;
    
  end %if nchan>2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % histogram_shift
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'histogram_shift')
  
  tmpcfg             = [];
  tmpcfg.method      = 'histogram';
  tmpcfg.equalbinavg = 'no';
  
  nshift = size(input{1},1);
  ncond  = length(input);
  %if nchan~=2, ft_error('number of trials ~= 2, do not know how to stratify'); end
  
  for k = 1:cfg.niter
    tmpinput = input;
    for cndlop = 1:ncond
      sel{cndlop}      = floor(abs(rand(size(input{cndlop},2),1)-eps)*nshift)+1;
      tmpinput{cndlop} = input{cndlop}(sel{cndlop}+6.*[0:length(sel{cndlop})-1]')';
    end
    
    [tmpoutput, binaxis] = stratify(tmpcfg, tmpinput{:});
    if k == 1
      bestoutput = tmpoutput;
      bestaxis   = binaxis;
      bestok     = sum(~isnan(tmpoutput{1}));
    end
    
    numok = sum(~isnan(tmpoutput{1}));
    [bestok, ind] = max([bestok numok]);
    if ind>1
      bestoutput = tmpoutput;
      bestaxis   = binaxis;
    end
  end
  
  %------create output
  varargout{1} = bestoutput;
  varargout{2} = bestaxis;
  varargout{3} = sel;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % splitxxxx
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(strfind(cfg.method, 'split'))
  ncond = length(varargin);
  if ~size(varargin{1},1)==2, ft_error('two channels required'); end
  m1      = mean(varargin{1},2);
  m2      = mean(varargin{2},2);
  sel{1} = zeros(1,size(varargin{1},2));
  sel{2} = zeros(1,size(varargin{2},2));
  
  sel{1}(find(varargin{1}(1,:) <= m1(1) & varargin{1}(2,:) <= m1(2))) = 1;
  sel{1}(find(varargin{1}(1,:) <= m1(1) & varargin{1}(2,:) >  m1(2))) = 2;
  sel{1}(find(varargin{1}(1,:) >  m1(1) & varargin{1}(2,:) <= m1(2))) = 3;
  sel{1}(find(varargin{1}(1,:) >  m1(1) & varargin{1}(2,:) >  m1(2))) = 4;
  
  sel{2}(find(varargin{2}(1,:) <= m2(1) & varargin{2}(2,:) <= m2(2))) = 1;
  sel{2}(find(varargin{2}(1,:) <= m2(1) & varargin{2}(2,:) >  m2(2))) = 2;
  sel{2}(find(varargin{2}(1,:) >  m2(1) & varargin{2}(2,:) <= m2(2))) = 3;
  sel{2}(find(varargin{2}(1,:) >  m2(1) & varargin{2}(2,:) >  m2(2))) = 4;
  
  if ~isempty(strfind(cfg.method, 'hilo'))
    sel{1} = sel{1} == 4;
    sel{2} = sel{2} == 1;
  elseif ~isempty(strfind(cfg.method, 'lohi'))
    sel{1} = sel{1} == 1;
    sel{2} = sel{2} == 4;
  elseif ~isempty(strfind(cfg.method, 'lolo'))
    sel{1} = sel{1} == 1;
    sel{2} = sel{2} == 1;
  elseif ~isempty(strfind(cfg.method, 'hihi'))
    sel{1} = sel{1} == 4;
    sel{2} = sel{2} == 4;
  end
  
  %------create output
  output = input;
  for cndlop = 1:ncond
    output{cndlop}(:,find(sel{cndlop}==0))=nan;
  end
  varargout{1} = output;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % equatespike
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'equatespike')
  if length(input)~=2
    ft_error('this requires two conditions as input');
  end
  datA = input{1};
  datB = input{2};
  
  ntrlA  = length(datA);
  ntrlB  = length(datB);
  nchanA = size(datA{1},1);
  nchanB = size(datB{1},1);
  nsmpA  = size(datA{1},2);
  nsmpB  = size(datB{1},2);
  
  if (ntrlA~=ntrlB)
    ft_error('number of trials should be the same in each condition');
  end
  
  if (nchanA~=nchanB)
    ft_error('number of channels should be the same in each condition');
  end
  
  if (nsmpA~=nsmpB)
    ft_error('number of samples should be the same in each condition');
  end
  
  for i=1:ntrlA
    if size(datA{i},2)~=nsmpA
      ft_error('number of samples should be the same in each trial');
    end
  end
  
  for i=1:ntrlB
    if size(datB{i},2)~=nsmpB
      ft_error('number of samples should be the same in each trial');
    end
  end
  
  nchan = nchanA;   % the same for each condition
  ntrl = ntrlA;     % the same for each condition
  
  if isequal(cfg.channel, 'all')
    chansel = 1:nchan;
  else
    chansel = cfg.channel;
  end
  
  % count the number of spikes in each trial and channel
  numA = zeros(ntrlA, nchanA);
  numB = zeros(ntrlB, nchanB);
  for i=1:ntrlA
    numA(i,:) = sum(datA{i},2);
  end
  for i=1:ntrlB
    numB(i,:) = sum(datB{i},2);
  end
  numA(:,setdiff(1:nchan, chansel)) = 0;
  numB(:,setdiff(1:nchan, chansel)) = 0;
  
  fprintf('average (over trials) of the number of spikes per channel in condition A: ');
  fprintf('%.0f ', mean(numA,1));
  fprintf('\n');
  fprintf('average (over trials) of the number of spikes per channel in condition B: ');
  fprintf('%.0f ', mean(numB,1));
  fprintf('\n');
  fprintf('standard deviation (over trials) of the number of spikes per channel in condition A: ');
  fprintf('%.0f ', std(numA,1));
  fprintf('\n');
  fprintf('standard deviation (over trials) of the number of spikes per channel in condition B: ');
  fprintf('%.0f ', std(numB,1));
  fprintf('\n');
  
  if strcmp(cfg.pairtrials, 'linkage')
    % concatenate the spike numbers for both conditions
    num = [numA; numB];
    % compute the distance between all possible pairs
    y = pdist(num, 'cityblock');
    
    % remove the connections between all pairs within a single condition
    y = squareform(y);
    sel = 1:ntrlA;
    y(sel,sel) = inf;
    sel = ntrlA + (1:ntrlB);
    y(sel,sel) = inf;
    for i=1:(ntrlA+ntrlB)
      y(i,i) = 0;
    end
    y = squareform(y);
    
    % determine the pairs that are the closest
    z = linkage(y, cfg.linkage);
    
    if any(any(z(1:ntrl,:)>2*ntrl))
      ft_error('trial pairs are not correct after hierarchical clustering');
    else
      fprintf('remaining distance after hierarchical clustering is %d\n', sum(z(1:ntrl,3)));
    end
    
    % ensure that the order of the pairs is always correct
    z = z(1:ntrl,:);
    for i=1:ntrl
      if z(i,1)>z(i,2)
        z(i,1:2) = z(i,2:1);
      end
    end
    
    % figure
    % hold on
    % for i=1:ntrl
    %   x = [num(z(i,1),1) num(z(i,2),1)];
    %   y = [num(z(i,1),2) num(z(i,2),2)];
    %   plot(x, y, '.-');
    % end
    
    indxA = z(:,1);
    indxB = z(:,2) - ntrlA;
    dist  = z(:,3);
    
  elseif strcmp(cfg.pairtrials, 'spikesort')
    [srtA, srtB, indxA, indxB] = spikesort(numA, numB, 'presort', 'global');
    dist = sum(abs(numA(indxA,:)-numB(indxB,:)),2);
    
  elseif strcmp(cfg.pairtrials, 'no')
    % no sorting to pair the trials is required
    indxA = 1:ntrlA;
    indxB = 1:ntrlA;
    dist = sum(abs(numA(indxA,:)-numB(indxB,:)),2);
    
  else
    ft_error('incorrect value for cfg.pairtrials');
  end % if pairtrials
  
  fprintf('removing %d spikes from a total of %d spikes\n', sum(dist), sum(sum(numA))+sum(sum(numB)));
  
  for i=1:ntrl
    delta = numA(indxA(i),:) - numB(indxB(i),:);
    outA{i} = datA{indxA(i)};
    outB{i} = datB{indxB(i)};
    for j=chansel(:)'
      if delta(j)>0
        % fprintf('removing %d spikes in condition A\n', abs(delta(j)));
        % make a vector (=rem) containing the indices of each spike. note that there can be multiple spikes at the same sample
        sel = find(outA{i}(j,:));
        cnt = outA{i}(j,sel);
        rem = zeros(sum(cnt), 1);
        s   = 0;
        for k=1:length(sel)
          for l=1:cnt(k)
            s = s + 1;
            rem(s) = sel(k);
          end
        end
        rem = rem(randperm(length(rem)));             % randomly permute the vector
        rem = rem(1:abs(delta(j)));                   % select random spikes to be removed
        for k=1:length(rem)
          outA{i}(j,rem(k)) = outA{i}(j,rem(k)) - 1;  % remove one spike from this sample
        end
      elseif delta(j)<0
        % fprintf('removing %d spikes in condition B\n', abs(delta(j)));
        % make a vector (=rem) containing the indices of each spike. note that there can be multiple spikes at the same sample
        sel = find(outB{i}(j,:));
        cnt = outB{i}(j,sel);
        rem = zeros(sum(cnt), 1);
        s   = 0;
        for k=1:length(sel)
          for l=1:cnt(k)
            s = s + 1;
            rem(s) = sel(k);
          end
        end
        rem = rem(randperm(length(rem)));             % randomly permute the vector
        rem = rem(1:abs(delta(j)));                   % select random spikes to be removed
        for k=1:length(rem)
          outB{i}(j,rem(k)) = outB{i}(j,rem(k)) - 1;  % remove one spike from this sample
        end
      else
        % spike numbers are equal, nothing to do
      end % if
    end % for nchan
  end % for ntrl
  varargout{1} = {outA, outB};
  
elseif strcmp(cfg.method, 'lohi')
  %%%%%%%%%%%%%%%%%
  %experimental code working on single channel inputs (2) selecting the
  %lowest amplitude for input 1 and highest for input 2
  %%%%%%%%%%%%%%%%%
  if length(varargin)~=2
    ft_error('two input arguments with data required');
  end
  if size(varargin{1},1)~=1 || size(varargin{2},1)~=1
    ft_error('only one channel per input is allowed');
  end
  if size(varargin{1},2)~=size(varargin{2},2)
    ft_error('the number of observations should be equal');
  end
  
  [srt1, ix1] = sort(input{1}, 'ascend');
  [srt2, ix2] = sort(input{2}, 'descend');
  mx1 = -inf;
  mx2 = inf;
  cnt = 0;
  sel = [];
  while mx2>mx1 && cnt<numel(srt1)
    cnt = cnt+1;
    sel = unique([sel ix1(cnt) ix2(cnt)]);
    mx1 = mean(input{1}(sel));
    mx2 = mean(input{2}(sel));
  end
  
  varargout{1} = input{1};
  varargout{2} = input{2};
  
  sel = setdiff(1:numel(srt1), sel);
  varargout{1}(sel) = nan;
  varargout{2}(sel) = nan;
  
end % cfg.method
