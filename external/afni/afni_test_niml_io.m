% a cell with filenames for testing
% for now these are all niml ascii files
fns=afni_niml_gettestfiles(); 

n=numel(fns);

for k=1:n %number 3 did not work as it's an ROI dataset
    fn=fns{k};
    
    [p,s]=afni_niml_read(fn);

    % print a few lines
    fprintf('%s\n', s(1:1000));
    
    % switch a few times between parsing and printing, results should
    % converge
    s2=afni_niml_print(p);
    p2=afni_niml_parse(s2);
    s3=afni_niml_print(p2);
    
    if ~isequal(s2,s3)
        error('Wrong string representation for %s', fn);
    end
    
    if ~isequal(p,p2)
        error('Wrong struct representation for %s', fn);
    end
    
    fprintf('Good, data seems to match\n');
    
    % make a new id code to keep suma happy
    p.self_idcode=['XYZ_' char(rand(1,24)*26+65)];
    
    % generate some random data, if there is any useful data (i.e. floats
    % as we don't want to mess with node indices)
    if isfield(p,'nodes')
        for j=1:numel(p.nodes)
            nk=p.nodes{j};
            if ~isempty(strfind(nk.ni_type,'float'))
                fprintf('Generating some random data\n');
                
                % ensure that we're in about the same range as the original
                % dataset (with some overshoot because of gaussianness)
                left=min(nk.data(:));
                right=max(nk.data(:));
                p.nodes{j}.data=left+(right-left)*randn(size(nk.data));
            end
        end
    end
        
    % convert to string
    [pth,f,e]=fileparts(fn);
    newfn=[pth '/niml_io_test_' f  e];
   
    afni_niml_write(newfn,p);
    
end
    