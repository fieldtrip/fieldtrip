function spm_print(job)
% Print the graphics window
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

% Run spm_print always as job to get configured print options
if nargin == 0
    spm_jobman('serial','','spm.util.print','');
    return;
elseif ischar(job)
    spm_jobman('serial','','spm.util.print',job);
    return;
end


try
    mon = {'Jan','Feb','Mar','Apr','May','Jun',...
            'Jul','Aug','Sep','Oct','Nov','Dec'};
    t   = clock;
    nam = ['spm_' num2str(t(1)) mon{t(2)} sprintf('%.2d',t(3))];

    if isempty(job.fname)
        if job.opts.append
            nam1 = fullfile(pwd,[nam job.opts.ext]);
        else
            nam1 = sprintf('%s_%3d',nam,1);
            for i=1:100000
                nam1 = fullfile(pwd,sprintf('%s_%.3d%s',nam,i,job.opts.ext));
                if ~exist(nam1,'file'), break; end;
            end
        end
    else
        nam1 = job.fname;
    end
    opts = {nam1,'-noui','-painters',job.opts.opt{:}};
    if strcmp(get(gcf,'Tag'),'Help'),
        fg = gcf;
    else
        fg = spm_figure('FindWin','Graphics');
    end
    if isdeployed
        deployprint(fg,opts{:});
    else
        print(fg,opts{:});
    end
    if isempty(strfind(nam1,filesep))
    fprintf('\nPrinting Graphics Windows to\n%s%s%s\n',pwd,filesep,nam1);
    else
    fprintf('\nPrinting Graphics Windows to\n%s\n',nam1);
    end
catch
    errstr = lasterror;
    errstr = errstr.message;
    str    = textscan(errstr,'%s','delimiter',sprintf('\n'));
    str    = str{1};
    str    = {str{:},'','- Print options are:', opts{:},...
                     '','- Current directory is:',['    ',pwd],...
                     '','            * nothing has been printed *'};
    spm('alert!',str,'printing problem...',sqrt(-1));
end
