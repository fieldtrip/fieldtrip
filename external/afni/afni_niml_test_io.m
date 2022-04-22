function afni_niml_test_io(fns)
% simple test function for i/o capabilities of afni_niml_{read,parse,print,write}
%
% AFNI_NIML_TEST_IO(FNS), where FNS is a cell with filenames, reads each of
% the files indicated in FNS, tests whether parsing and printing are
% reverse operations, and writes the files to a new file. If a file
% contains a matrix with floats, then these are replaced by random numbers.
%
% AFNI_NIML_TEST_IO(FN), where FN is a single filename, is also allowed.
% AFNI_NIML_TEST_IO() uses default filenames as defined in the body of this
% function.
%
% Each file named in FNS should be in ASCII NIML format.
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

if nargin<1
    % set defaults
    fns={'testfiles/you_look_marvellous.niml.dset',...
        'testfiles/test_ROI2.niml.dset',...
        'testfiles/moredata.niml.dset'};
end

if ischar(fns) % single filename
    fn=fns;
    fns=cell(1);
    fns{1}=fn;
end

n=numel(fns);

for k=1:n
    fn=fns{k};
    fprintf('\n*** Processing %s ***\n\n',fn);



    [p,s]=afni_niml_read(fn);

    % print a few lines
    fprintf('First 1000 characters:\n%s\n', s(1:1000));

    % switch a few times between parsing and printing, results should
    % converge to 'fixed points'
    fprintf('Check whether parsing and printing are mutually inverse\n');
    s2=afni_niml_print(p);
    p2=afni_niml_parse(s2);
    s3=afni_niml_print(p2);

    % check convergence
    if ~isequal(s2,s3)
        error('Wrong string representation for %s', fn);
    end

    if ~isequal(p,p2)
        error('Wrong struct representation for %s', fn);
    end

    fprintf('Good, data seems to match\n');

    % now read data as a 'simple' struct, that is one that we can
    % manipulate easily
    S=afni_niml_readsimple(fn);

    fprintf('Simple struct read:\n');
    disp(S);

    % print the information for each column (a la AFNI SurfDsetinfo)
    [nverts,ncols]=size(S.data);
    fprintf('Found %d columns\n', ncols);
    for j=1:ncols
        fprintf('COlumn %d: label "%s" with stat value "%s"\n', j, S.labels{j}, S.stats{j});
    end

    % for every odd colum, make values gaussian random
    fprintf('Put some random data in\n');
    for j=1:2:ncols
        d=S.data(:,j);
        if isequal(d,round(d)) % if integers, ignore (probably ROI indices)
            continue;
        end
        S.data(:,j)=randn(nverts,1);
        S.labels{j}=sprintf('rnd(%d)',j);
        S.stats{j}='Zscore()';
    end

    % convert to string
    [pth,f e]=fileparts(fn);
    newfn=[pth '/niml_io_test_' f e];

    afni_niml_writesimple(newfn,S);
end


% Make a file from scratch
myrootfn='testfiles/myfile';
nverts=100002;
S=struct();
S.labels={'Hello','World'};
S.stats={'Zscore()','none'};
S.data=[randn(nverts,1) rand(nverts,1)];

% write with all info
afni_niml_writesimple(S,[myrootfn '_1.niml.dset']);

% just write the data
afni_niml_writesimple(S.data,[myrootfn '_2.niml.dset']);


