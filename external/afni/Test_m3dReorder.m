%set the filename identifier
Identifiers = {'AAzst1avir.N02.*.HEAD', 'AAzst1avir.N05.*.HEAD'}; %Modify here
[err, ErrMessage, List] = zglobb (Identifiers);

%loop across all Bricks found
Nel = length(List);
for (i=1:1:Nel),
	Input = List(i).name; 
	fprintf(1,'\nNow processing: %s ...', Input);
	[ans, I_Prefix, View] = PrefixStatus(Input);
	%set the new prefix
	Prefix = sprintf('%s_reord', I_Prefix); %Modify here
	%set up for the function m3dReorder
	Mapfile = 'map.1D'; %Modify here
	Opt.Verbose = 1; %Modify here
	Opt.Detrend = 2; %Modify here
	Opt.Dup = 'Col'; %Modify here
	Opt.NoCheck = 0;
	[err] = m3dReorder (Input, Prefix, Mapfile, Opt);
end
