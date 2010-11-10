function ExamineXmat(fname, polort, dt, nrun, nmot)
%  ExamineXmat([XMAT])
%
%  A function to aid in examining the design matrix output by 3dDeconvolve
%  
% XMAT: Name of xmat file, such as EA.xmat.1D. If unspecified,
%        function will open a file browser to pick one.
% The function will begin by displaying all the columns (regressors)
% in XMAT.
% The title of the graph shows:
%   name of the file
%   TR
%   condition numbers for:
%      entire matrix (Full)
%      entire matrix without motion regressors (NoMot)
%      entire matrix without motion and without baseline (NoMotNoBase)
%      plotted columns only (Viewed)
%
% Labels of various columns are printed to the left, if there are not
% too many columns displayed. Labels contain the column index and the label, 
% separated by a ":"
%
% The function will prompt you for column selections, which would allow
% you to examine a section of the design matrix at a time. Looking at 
% the 'Viewed' condition number would help you find which task
% regressors may be causing multicollinearity. Column selection 
% can be done using strings that match one or a set of regressors or 
% using column indices. The prompt interface has selection examples. 
%

% old use: ExamineXmat(fname, polort, dt, nrun, nmot)
verb = 0;

if (nargin < 1 | isempty(fname) | (ischar(fname) & ~filexist(fname))),
   [fname, pname] = uigetfile('*.xmat.1D','Pick an Xmat');
   fname = sprintf('%s%c%s', pname, filesep, fname);
   if (~filexist(fname)),
      fprintf(2,'File %s not found\n', fname);
      return;
   end
end
if (nargin < 2),
   csstims = cellstr('');
   polort = -1;
   nrun = 0;
   cntstims = 0;
   %get the info from 3dSynthesize
   if (verb) fprintf(2,'Running 3dSynthesize...\n'); end
   com = sprintf('3dSynthesize -dry_info -matrix %s -select all', fname);
   [s1, s2] = unix(com);
   if (verb) fprintf(2,'Parsing...\n'); end
   cnt = 1;
   imot = [];
   while (~isempty(s2)),
      [d,c,e,n] = sscanf(s2,'%s ',1);
      s2 = s2(n:length(s2));
      if (strncmp(d,'TR:',3)),
         [tt,dt] = strread(d,'%s%f','delimiter',':');
      elseif (~isempty(strfind(d,'Run#'))),
         [col,tt] = strread(d,'%d%s','delimiter',':'); 
         if (strncmp(tt,'Run#1',5)), 
            polort = polort+1;   
         elseif (strncmp(tt,'Run#',4)),
            nrun = nrun+1;
         elseif (strncmp(tt(length(tt):-1:1),'0#',2)),
            cntstims = cntstims+1;
            csstims(cntstims) = cellstr(tt);
         end 
         cs(cnt) = cellstr(d);
         cnt = cnt + 1;
      elseif ( ~isempty(strfind(d,'roll#'))  ||...
               ~isempty(strfind(d,'pitch#')) ||...
               ~isempty(strfind(d,'yaw#'))   ||...
               ~isempty(strfind(d,'dS#'))    ||...
               ~isempty(strfind(d,'dL#'))    ||...
               ~isempty(strfind(d,'dP#'))    ||...
               ~isempty(strfind(d,'x-scale#'))||...
               ~isempty(strfind(d,'y-scale#'))||...
               ~isempty(strfind(d,'z-scale#'))||...
               ~isempty(strfind(d,'yx-shear#'))||...
               ~isempty(strfind(d,'zx-shear#'))||...
               ~isempty(strfind(d,'zy-shear#'))),
         [col,tt] = strread(d,'%d%s','delimiter',':'); 
         imot = [imot col];
         cs(cnt) = cellstr(d);
         cnt = cnt + 1;
      elseif (~isempty(strfind(d,'#'))),
         if (strncmp(d([length(d):-1:1]),'0#',2)),
            cntstims = cntstims+1;
            [col,tt] = strread(d,'%d%s','delimiter',':'); 
            csstims(cntstims) = cellstr(tt);
         end 
         cs(cnt) = cellstr(d);
         cnt = cnt + 1;
      end
   end
   nrun = nrun/(polort+1)+1;
   nmot = length(imot);
end
if (polort < 0),
   polort = [];
   while(isempty(polort) | ~isnumeric(polort) | polort < 0),
      polort = input('Enter polort: ');
   end
end
if (dt < 0.0),   
   dt = [];
   while(isempty(dt) | ~isnumeric(dt) | dt < 0.0),
      dt = input('Enter dt: ');
   end
end
if (nrun < 0.0),   
   nrun = [];
   while(isempty(nrun) | ~isnumeric(nrun) | nrun < 0.0),
      nrun = input('Enter number of runs: ');
   end
end
if (nmot < 0.0),   
   nmot = [];
   while(isempty(nmot) | ~isnumeric(nmot) | nmot < 0.0),
      nmot = input('Enter number of motion regressors (assumed at the very end): ');
   end
end

if (verb) fprintf(2,'Reading 1D...\n'); end
Opt.verb=verb;
Opt.method = 3;
if (ischar(fname)) [e,Xabi] = Read_1D(fname, Opt);
else Xabi = fname; fname = 'matrix in mem.';
end

s = '---';

%remove baseline and, assuming motion is last, remove last 6 regressors
%but we do not reliably know where motion lies, so for now, default is nmot = 0

istim = [1+(polort+1)*nrun, size(Xabi,2)-nmot];
trimmed = Xabi(:,[istim(1):1:istim(2)]);
CondFull = cond(Xabi);
CondNoMotion = cond(Xabi(:, [1:size(Xabi,2)-nmot]));
CondNoMotionNoBase = cond(Xabi(:,[1+(polort+1)*nrun:size(Xabi,2)-nmot]));

mshow = Xabi;  
UseMult = 1;
dogui = 1;

if (verb) fprintf(2,'Display time...\n'); end

v = [1:1:size(mshow,2)];
t = 0:dt:(size(mshow,1)-1)*dt;
while (~isempty(s)),
   if (isempty(v)), %maybe bad parsing
      snote = sprintf('No proper selection, showing entire matrix.');
      if (dogui),
         warndlg(snote);
      else
         fprintf(2,'\n%s\n', snote);
      end
      v = [1:1:size(Xabi,2)];
   end
   inz = [];
   deltaoffs = 2;
   figure(1); clf
   if (length(v) < 10 & UseMult),
      for (i=1:1:length(v)),
         if (i==1),
            titplot = subplot (length(v), 1, i);
         else
            subplot (length(v), 1, i);
         end
         if (length(v) < 7),
            plot (t, mshow(:, v(i)), 'b-', t, mshow(:, v(i)), 'r*');
            if (sum(abs(mshow(:,v(i)))) == 0),
               plot (t(1:max([1,length(t)/20]):end),...
                      mshow(1:max([1,length(t)/20]):end,v(i))+offs,...
                       'bo'); hold on
            end
            ylabel(char(cs(v(i))), 'interpreter', 'none');
         else
            plot (t, mshow(:, v(i)), 'b-');
            if (sum(abs(mshow(:,v(i)))) == 0),
               plot (t(1:max([1,length(t)/20]):end),...
                      mshow(1:max([1,length(t)/20]):end,v(i))+offs,...
                       'bo'); hold on
            else
               inz = [inz v(i)];
            end
            ylabel(char(cs(v(i))), 'interpreter', 'none');
         end
         if (~isempty(cs) & length(v)<20),
            %ylabel(char(cs(v(i))), 'interpreter', 'none');
            str = char(cs(v(i)));
         else 
            %ylabel(sprintf('R%d', v(i)+1),  'interpreter', 'none');
            str = sprintf('R%d', v(i)+1);
         end
         ms = mshow(:,v(i));
         mpp = ms(find(ms > 0));
         if (~isempty(mpp)),
            deltaoffs = mean(mpp);
            offs = offs+deltaoffs;
         else
            offs = offs+deltaoffs;
         end
         if (length(v) < 35),
            xt = get(gca,'XTick');
            text( xt(1)-0.1.*xt(2), offs, char(cs(v(i))), ...
                  'HorizontalAlignment', 'right',...
                  'interpreter', 'none');
         end
      end
   else
      offs = 0;
      titplot = subplot (1,1,1);
      for (i=1:1:length(v)),
         ms = mshow(:,v(i));
         mpp = ms(find(ms > 0));
         if (~isempty(mpp)),
            deltaoffs = mean(mpp);
            offs = offs+deltaoffs;
         else
            offs = offs+deltaoffs;
         end
         ccc = ROIcol;
         plot (t, mshow(:,v(i))+offs, 'Color', ccc); hold on
         if (sum(abs(mshow(:,v(i)))) == 0),
            plot (t(1:max([1,length(t)/20]):end),...
                   mshow(1:max([1,length(t)/20]):end,v(i))+offs,...
                    'o', 'Color', ccc); hold on
         else
            inz = [inz v(i)];
         end
         if (length(v) < 75),
            xt = get(gca,'XTick');
            text( xt(1)-0.1.*xt(2), offs, char(cs(v(i))), ...
                  'HorizontalAlignment', 'right',...
                  'interpreter', 'none');
         end
      end
      set(gca,'YTickLabel',[]);
   end
   xlabel('time (sec)');
   subplot(titplot);
   title(sprintf(['',...
            'Xmat %s,   TR %.3f,    %d regressors in view,   query >> %s <<\n'...
            'Condition #: '...
            'Full %d\t NoMot %g\t NoMotNoBase %g \t'...
            'Viewed %g, Viewed no zeros %g\n'],...
             fname, dt, length(v), s, CondFull, CondNoMotion,...
             CondNoMotionNoBase, cond(mshow(:, [v])),...
             cond(mshow(:, [inz]))),...
            'interpreter', 'none'); 
   
   if (1), %show correlation matrix for Assaf
      figure(2); clf
      ccm = corrcoef(mshow(:,v)); 
      ccm = ccm.*(1-diag(diag(ccm,0)));
      mmm = max(abs(ccm(:)));
      imagesc(ccm, [-mmm , mmm]); colorbar;
      ima = gca; 
      title('Correlation matrix of selected regressors. Diagonal set to 0');
      set (ima, 'YTick', [1:1:size(ccm,1)]');
      set (ima, 'YTickLabel', cellstr(cs(v)));
      set (ima, 'XTick', [1:1:size(ccm,1)]');
      set (ima, 'XTickLabel', cellstr(cs(v)));
      xticklabel_rotate([], 90);
   end
   
   smpl = char(csstims(1));
   if (length(csstims)>2),
      sss =  char(csstims(length(csstims)));
      smpl2 = sprintf(' or "%s %s"',...
                      smpl(1:min(4,length(smpl))),...
                      sss(1:min(4,length(sss))));
   else
      smpl2 = '';
   end
   if (istim(2) - istim(1) > 5),
      smplcol = sprintf('"%d", "[%d,%d,%d]", or "[%d:%d]"',...
                        istim(1)-1,...
                        istim(1), istim(1)+2, istim(2)-1,...
                        istim(1)-1, istim(2)-1);
   else,
      smplcol = sprintf('"3", "[0,1]", or "[2:4]"');
   end
   sprompt = sprintf([ 'Enter which regressors (columns) you want to see. \n'...
                       '-------------------------------------------------\n'...
                       'You can use text such as "%s"%s that matches one or '...
                       'a few labels,\n'...
                       'or column indices in matlab''s vector notation '...
                       'such as %s\n'...
                       '--\n',...
                       'Column indices start at 0 as with AFNI.\n'...
                       '--\n',...
                       'In this XMAT stimuli are in columns between %d and %d:\n '...
                       '--\n'...
                       ], ...
                       smpl(1:min(4,length(smpl))),...
                       smpl2,...
                       smplcol,...
                       istim(1)-1, istim(2)-1);
   bads = 1;
   while (bads),
      if (~dogui),
         s = zdeblank(input (sprompt, 's'));
      else,
         Opt.Resize='on';
         Opt.WindowStyle='normal';
         s = zdeblank(char(zinputdlg(cellstr(sprompt),...
                                    'Select Regressors',...
                                    1,...
                                    cellstr(''),...
                                    Opt)));
      end 
      ipunc = find (s == ',');
      s(ipunc) = ' ';
      if (isempty(s)), 
         return;
      else
         if (isdigit(s(1)) | s(1) == '[' | s(1) == '(')
            eval(sprintf('v=%s;', s)); 
            if (0 & isempty(v)),
               return;
            end
            v = v + 1;
         else
            v = [];
            stmp = s;
            while (~isempty(stmp)),
               [d,c,e,n] = sscanf(stmp,'%s ',1);
               stmp = stmp(n:length(stmp));
               d = zdeblank(d);
               for(i=1:1:length(cs)),
                  [col,tt] = strread(char(cs(i)),'%d%s','delimiter',':');
                  if (strncmp(tt,d, length(d))),
                     v = [v i];
                  end
               end
            end
            if (0 & isempty(v)),
               return;
            end
         end
      end
      if (length(v) && (min(v) < 1 || max(v) > size(Xabi,2))),
         if (dogui),
            warndlg(sprintf('Bad range for column selection:\nmin=%d, max=%d\nlimit is min=%d, max=%d\n', min(v), max(v), 1, size(Xabi,2)));
         else
            fprintf(2,'Bad range for column selection:\nmin=%d, max=%d\nlimit is min=%d, max=%d\n', min(v), max(v), 1, size(Xabi,2));
         end
         bads = 1;
      else
         bads = 0;
      end
   end
end

