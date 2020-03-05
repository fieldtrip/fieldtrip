function S = read_mclust_t(tfilelist)

% adapted from M-clust function LoadSpikes

%-------------------
% Check input type
%-------------------
if ~isa(tfilelist, 'cell')
   ft_error('LoadSpikes: tfilelist should be a cell-array.');
end

nFiles = length(tfilelist);


S = cell(nFiles, 1);
for iF = 1:nFiles
	tfn = tfilelist{iF};
	if ~isempty(tfn)
    try
 		  tfp = fopen_or_error(tfn, 'rb','b');
    catch err
			warning([ 'Could not open tfile ' tfn]);
      continue
		end
		
		ReadHeader(tfp);    
		S{iF} = fread(tfp,inf,'uint64');	
	  S{iF} = double(S{iF}*100);

		fclose(tfp);		
	end 		% if tfn valid
end		% for all files
fprintf(2,'\n');

