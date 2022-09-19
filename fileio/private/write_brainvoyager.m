function write_brainvoyager(filename, data, dataformat, vmpversion)

% helper function to write volumetric data for brainvoyager.
% this is old code that moved from ft_volumewrite to clean up
% the high level function a bit. it is assumed that the orientation
% of the volume is correct.

data = double(data);
siz  = size(data);
maxval = max(data(:));

[p, f, e] = fileparts(filename);
filename  = fullfile(p, f);

switch dataformat
  case 'vmp'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in BrainVoyager VMP format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(sprintf('%s.vmp', filename),'w');
    if fid < 0
      ft_error('Cannot write to file %s.vmp\n',filename);
    end
    
    switch vmpversion
      case 1
        % write the header
        fwrite(fid, 1, 'short');      % version
        fwrite(fid, 1, 'short');      % number of maps
        fwrite(fid, 1, 'short');      % map type
        fwrite(fid, 0, 'short');      % lag
        
        fwrite(fid, 0, 'short');      % cluster size
        fwrite(fid, 1, 'float');      % thresh min
        fwrite(fid, maxval, 'float'); % thresh max
        fwrite(fid, 0, 'short');      % df1
        fwrite(fid, 0, 'short');      % df2
        fwrite(fid, 0, 'char');       % name
        
        fwrite(fid, siz, 'short');    % size
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(1)-1, 'short');
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(2)-1, 'short');
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(3)-1, 'short');
        fwrite(fid, 1, 'short');      % resolution
        
        % write the data
        fwrite(fid, data, 'float');
      case 2
        % determine relevant subvolume
        % FIXME, this is not functional at the moment, since earlier in this function all nans have been replaced by zeros
        minx = min(find(~isnan(max(max(data,[],3),[],2))));
        maxx = max(find(~isnan(max(max(data,[],3),[],2))));
        miny = min(find(~isnan(max(max(data,[],3),[],1))));
        maxy = max(find(~isnan(max(max(data,[],3),[],1))));
        minz = min(find(~isnan(max(max(data,[],1),[],2))));
        maxz = max(find(~isnan(max(max(data,[],1),[],2))));
        
        % write the header
        fwrite(fid, 2, 'short');      % version
        fwrite(fid, 1, 'int');        % number of maps
        fwrite(fid, 1, 'int');        % map type
        fwrite(fid, 0, 'int');        % lag
        
        fwrite(fid, 0, 'int');        % cluster size
        fwrite(fid, 0, 'char');       % cluster enable
        fwrite(fid, 1, 'float');      % thresh
        fwrite(fid, maxval, 'float'); % thresh
        fwrite(fid, 0, 'int');        % df1
        fwrite(fid, 0, 'int');        % df2
        fwrite(fid, 0, 'int');        % bonf
        fwrite(fid, [255,0,0], 'uchar');   % col1
        fwrite(fid, [255,255,0], 'uchar'); % col2
        fwrite(fid, 1, 'char');       % enable SMP
        fwrite(fid, 1, 'float');      % transparency
        fwrite(fid, 0, 'char');       % name
        
        fwrite(fid, siz, 'int');      % original size
        fwrite(fid, minx-1, 'int');
        fwrite(fid, maxx-1, 'int');
        fwrite(fid, miny-1, 'int');
        fwrite(fid, maxy-1, 'int');
        fwrite(fid, minz-1, 'int');
        fwrite(fid, maxz-1, 'int');
        fwrite(fid, 1, 'int');        % resolution
        
        % write the data
        fwrite(fid, data(minx:maxx,miny:maxy,minz:maxz), 'float');
    end
    fclose(fid);
    
  case 'vmr'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in BrainVoyager VMR format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(sprintf('%s.vmr',filename),'w');
    if fid < 0
      ft_error('Cannot write to file %s.vmr\n',filename);
    end
    
    % data should be scaled between 0 and 225
    data = data - min(data(:));
    data = round(225*data./max(data(:)));
    
    % write the header
    fwrite(fid, siz, 'ushort');
    % write the data
    fwrite(fid, data, 'uint8');
    fclose(fid);
end