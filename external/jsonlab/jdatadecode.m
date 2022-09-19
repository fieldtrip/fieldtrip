function newdata=jdatadecode(data,varargin)
%
% newdata=jdatadecode(data,opt,...)
%
% Convert all JData object (in the form of a struct array) into an array
% (accepts JData objects loaded from either loadjson/loadubjson or 
% jsondecode for MATLAB R2016b or later)
%
% This function implements the JData Specification Draft 3 (Jun. 2020)
% see http://github.com/fangq/jdata for details
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%      data: a struct array. If data contains JData keywords in the first
%            level children, these fields are parsed and regrouped into a
%            data object (arrays, trees, graphs etc) based on JData 
%            specification. The JData keywords are
%               "_ArrayType_", "_ArraySize_", "_ArrayData_"
%               "_ArrayIsSparse_", "_ArrayIsComplex_", 
%               "_ArrayZipType_", "_ArrayZipSize", "_ArrayZipData_"
%      opt: (optional) a list of 'Param',value pairs for additional options 
%           The supported options include
%               Recursive: [1|0] if set to 1, will apply the conversion to 
%                            every child; 0 to disable
%               Base64: [0|1] if set to 1, _ArrayZipData_ is assumed to
%                         be encoded with base64 format and need to be
%                         decoded first. This is needed for JSON but not
%                         UBJSON data
%               Prefix: ['x0x5F'|'x'] for JData files loaded via loadjson/loadubjson, the
%                         default JData keyword prefix is 'x0x5F'; if the
%                         json file is loaded using matlab2018's
%                         jsondecode(), the prefix is 'x'; this function
%                         attempts to automatically determine the prefix;
%                         for octave, the default value is an empty string ''.
%               FullArrayShape: [0|1] if set to 1, converting _ArrayShape_ 
%                         objects to full matrices, otherwise, stay sparse
%               MaxLinkLevel: [0|int] When expanding _DataLink_ pointers,
%                         this sets the maximum level of recursion
%               FormatVersion: [2|float]: set the JSONLab output version; 
%                         since v2.0, JSONLab uses JData specification Draft 1
%                         for output format, it is incompatible with all
%                         previous releases; if old output is desired,
%                         please set FormatVersion to 1
%
% output:
%      newdata: the covnerted data if the input data does contain a JData 
%               structure; otherwise, the same as the input.
%
% examples:
%      obj={[],{'test'},true,struct('sparse',sparse(2,3),'magic',uint8(magic(5)))}
%      jdata=jdatadecode(jdataencode(obj))
%      isequaln(obj,jdata)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

    newdata=data;
    opt=struct;
    if(nargin==2)
        opt=varargin{1};
    elseif(nargin>2)
        opt=varargin2struct(varargin{:});
    end
    opt.fullarrayshape=jsonopt('FullArrayShape',0,opt);
    opt.maxlinklevel=jsonopt('MaxLinkLevel',0,opt);

    %% process non-structure inputs
    if(~isstruct(data))
        if(iscell(data))
            newdata=cellfun(@(x) jdatadecode(x,opt),data,'UniformOutput',false);
        elseif(isa(data,'containers.Map'))
            newdata=containers.Map('KeyType',data.KeyType,'ValueType','any');
            names=data.keys;
            for i=1:length(names)
                newdata(names{i})=jdatadecode(data(names{i}),opt);
            end
        end
        return;
    end
    
    %% assume the input is a struct below
    fn=fieldnames(data);
    len=length(data);
    needbase64=jsonopt('Base64',0,opt);
    format=jsonopt('FormatVersion',2,opt);
    if(isoctavemesh)
        prefix=jsonopt('Prefix','',opt);
    else
        prefix=jsonopt('Prefix','x0x5F',opt);
    end
    if(~isfield(data,N_('_ArrayType_')) && isfield(data,'x_ArrayType_'))
        prefix='x';
        opt.prefix='x';
    end

    %% recursively process subfields
    if(jsonopt('Recursive',1,opt)==1)
      for i=1:length(fn) % depth-first
        for j=1:len
            if(isstruct(data(j).(fn{i})) || isa(data(j).(fn{i}),'containers.Map'))
                newdata(j).(fn{i})=jdatadecode(data(j).(fn{i}),opt);
            elseif(iscell(data(j).(fn{i})))
                newdata(j).(fn{i})=cellfun(@(x) jdatadecode(x,opt),newdata(j).(fn{i}),'UniformOutput',false);
            end
        end
      end
    end

    %% handle array data
    if(isfield(data,N_('_ArrayType_')) && (isfield(data,N_('_ArrayData_')) || isfield(data,N_('_ArrayZipData_'))))
      newdata=cell(len,1);
      for j=1:len
        if(isfield(data,N_('_ArrayZipSize_')) && isfield(data,N_('_ArrayZipData_')))
            zipmethod='zip';
            if(isstruct(data(j).(N_('_ArrayZipSize_'))))
                data(j).(N_('_ArrayZipSize_'))=jdatadecode(data(j).(N_('_ArrayZipSize_')),opt);
            end
	    dims=data(j).(N_('_ArrayZipSize_'))(:)';
            if(length(dims)==1)
                 dims=[1 dims];
            end
            if(isfield(data,N_('_ArrayZipType_')))
                zipmethod=data(j).(N_('_ArrayZipType_'));
            end
            if(ismember(zipmethod,{'zlib','gzip','lzma','lzip','lz4','lz4hc','base64'}))
                decompfun=str2func([zipmethod 'decode']);
                arraytype=data(j).(N_('_ArrayType_'));
                chartype=0;
                if(strcmp(arraytype,'char') || strcmp(arraytype,'logical'))
                    chartype=1;
                    arraytype='uint8';
                end
                if(needbase64 && strcmp(zipmethod,'base64')==0)
                    ndata=reshape(typecast(decompfun(base64decode(data(j).(N_('_ArrayZipData_')))),arraytype),dims);
                else
                    ndata=reshape(typecast(decompfun(data(j).(N_('_ArrayZipData_'))),arraytype),dims);
                end
                if(chartype)
                    ndata=char(ndata);
                end
            else
                error('compression method is not supported');
            end
        else
            if(isstruct(data(j).(N_('_ArrayData_'))))
                data(j).(N_('_ArrayData_'))=jdatadecode(data(j).(N_('_ArrayData_')),opt);
            end
            if(isstruct(data(j).(N_('_ArrayData_'))) && isfield(data(j).(N_('_ArrayData_')),N_('_ArrayType_')))
                data(j).(N_('_ArrayData_'))=jdatadecode(data(j).(N_('_ArrayData_')),varargin{:});
            end
            if(iscell(data(j).(N_('_ArrayData_'))))
                data(j).(N_('_ArrayData_'))=cell2mat(cellfun(@(x) double(x(:)),data(j).(N_('_ArrayData_')),'uniformoutput',0)).';
            end
            ndata=cast(data(j).(N_('_ArrayData_')),char(data(j).(N_('_ArrayType_'))));
        end
        if(isfield(data,N_('_ArrayZipSize_')))
            if(isstruct(data(j).(N_('_ArrayZipSize_'))))
                data(j).(N_('_ArrayZipSize_'))=jdatadecode(data(j).(N_('_ArrayZipSize_')),opt);
            end
	    dims=data(j).(N_('_ArrayZipSize_'))(:)';
            if(iscell(dims))
                dims=cell2mat(dims);
            end
            if(length(dims)==1)
                 dims=[1 dims];
            end
            ndata=reshape(ndata(:),fliplr(dims));
            ndata=permute(ndata,ndims(ndata):-1:1);
        end
        iscpx=0;
        if(isfield(data,N_('_ArrayIsComplex_')) && isstruct(data(j).(N_('_ArrayIsComplex_'))) )
                data(j).(N_('_ArrayIsComplex_'))=jdatadecode(data(j).(N_('_ArrayIsComplex_')),opt);
        end
        if(isfield(data,N_('_ArrayIsComplex_')) && data(j).(N_('_ArrayIsComplex_')) )
                iscpx=1;
        end
        iscol=0;
        if(isfield(data,N_('_ArrayOrder_')))
            arrayorder=data(j).(N_('_ArrayOrder_'));
            if(~isempty(arrayorder) && (arrayorder(1)=='c' || arrayorder(1)=='C'))
                iscol=1;
            end
        end
        if(isfield(data,N_('_ArrayIsSparse_')) && isstruct(data(j).(N_('_ArrayIsSparse_'))) )
                data(j).(N_('_ArrayIsSparse_'))=jdatadecode(data(j).(N_('_ArrayIsSparse_')),opt);
        end
        if(isfield(data,N_('_ArrayIsSparse_')) && data(j).(N_('_ArrayIsSparse_')))
                if(isfield(data,N_('_ArraySize_')))
                    if(isstruct(data(j).(N_('_ArraySize_'))))
                        data(j).(N_('_ArraySize_'))=jdatadecode(data(j).(N_('_ArraySize_')),opt);
                    end
                    dim=data(j).(N_('_ArraySize_'))(:)';
                    if(iscell(dim))
                        dim=cell2mat(dim);
                    end
                    dim=double(dim);
		    if(length(dim)==1)
			dim=[1 dim];
		    end
                    if(iscpx)
                        ndata(end-1,:)=complex(ndata(end-1,:),ndata(end,:));
                    end
                    if isempty(ndata)
                        % All-zeros sparse
                        ndata=sparse(dim(1),prod(dim(2:end)));
                    elseif dim(1)==1
                        % Sparse row vector
                        ndata=sparse(1,ndata(1,:),ndata(2,:),dim(1),prod(dim(2:end)));
                    elseif dim(2)==1
                        % Sparse column vector
                        ndata=sparse(ndata(1,:),1,ndata(2,:),dim(1),prod(dim(2:end)));
                    else
                        % Generic sparse array.
                        ndata=sparse(ndata(1,:),ndata(2,:),ndata(3,:),dim(1),prod(dim(2:end)));
                    end
                else
                    if(iscpx && size(ndata,2)==4)
                        ndata(3,:)=complex(ndata(3,:),ndata(4,:));
                    end
                    ndata=sparse(ndata(1,:),ndata(2,:),ndata(3,:));
                end
        elseif(isfield(data,N_('_ArrayShape_')))
                if(isstruct(data(j).(N_('_ArrayShape_'))))
                    data(j).(N_('_ArrayShape_'))=jdatadecode(data(j).(N_('_ArrayShape_')),opt);
                end
                if(iscpx)
                    if(size(ndata,1)==2)
                        dim=size(ndata);
                        dim(end+1)=1;
                        arraydata=reshape(complex(ndata(1,:),ndata(2,:)),dim(2:end));
                    else
                        error('The first dimension must be 2 for complex-valued arrays');
                    end
                else
                    arraydata=data.(N_('_ArrayData_'));
                end
                shapeid=data.(N_('_ArrayShape_'));
                if(isfield(data,N_('_ArrayZipSize_')))
                        datasize=data.(N_('_ArrayZipSize_'));
                        if(iscell(datasize))
                            datasize=cell2mat(datasize);
                        end
                        datasize=double(datasize);
                        if(iscpx)
                            datasize=datasize(2:end);
                        end
                else
                        datasize=size(arraydata);
                end
                if(isstruct(data(j).(N_('_ArraySize_'))))
                    data(j).(N_('_ArraySize_'))=jdatadecode(data(j).(N_('_ArraySize_')),opt);
                end
                arraysize=data.(N_('_ArraySize_'));

                if(iscell(arraysize))
                    arraysize=cell2mat(arraysize);
                end
                arraysize=double(arraysize);
                if(ischar(shapeid))
                        shapeid={shapeid};
                end
                arraydata=double(arraydata).';
                if(strcmpi(shapeid{1},'diag'))
                        ndata=spdiags(arraydata(:),0,arraysize(1),arraysize(2));
                elseif(strcmpi(shapeid{1},'upper') || strcmpi(shapeid{1},'uppersymm'))
                        ndata=zeros(arraysize);
                        ndata(triu(true(size(ndata)))')=arraydata(:);
                        if(strcmpi(shapeid{1},'uppersymm'))
                            ndata(triu(true(size(ndata))))=arraydata(:);
                        end
                        ndata=ndata.';
                elseif(strcmpi(shapeid{1},'lower') || strcmpi(shapeid{1},'lowersymm'))
                        ndata=zeros(arraysize);
                        ndata(tril(true(size(ndata)))')=arraydata(:);
                        if(strcmpi(shapeid{1},'lowersymm'))
                            ndata(tril(true(size(ndata))))=arraydata(:);
                        end
                        ndata=ndata.';
                elseif(strcmpi(shapeid{1},'upperband') || strcmpi(shapeid{1},'uppersymmband'))
                        if(length(shapeid)>1 && isvector(arraydata))
                            datasize=double([shapeid{2}+1, prod(datasize)/(shapeid{2}+1)]);
                        end
                        ndata=spdiags(reshape(arraydata,min(arraysize),datasize(1)),-datasize(1)+1:0,arraysize(2),arraysize(1)).';
                        if(strcmpi(shapeid{1},'uppersymmband'))
                            diagonal=diag(ndata);
                            ndata=ndata+ndata.';
                            ndata(1:arraysize(1)+1:end)=diagonal;
                        end
                elseif(strcmpi(shapeid{1},'lowerband') || strcmpi(shapeid{1},'lowersymmband'))
                        if(length(shapeid)>1 && isvector(arraydata))
                            datasize=double([shapeid{2}+1, prod(datasize)/(shapeid{2}+1)]);
                        end
                        ndata=spdiags(reshape(arraydata,min(arraysize),datasize(1)),0:datasize(1)-1,arraysize(2),arraysize(1)).';
                        if(strcmpi(shapeid{1},'lowersymmband'))
                            diagonal=diag(ndata);
                            ndata=ndata+ndata.';
                            ndata(1:arraysize(1)+1:end)=diagonal;
                        end
                elseif(strcmpi(shapeid{1},'band'))
                        if(length(shapeid)>1 && isvector(arraydata))
                            datasize=double([shapeid{2}+shapeid{3}+1, prod(datasize)/(shapeid{2}+shapeid{3}+1)]);
                        end
                        ndata=spdiags(reshape(arraydata,min(arraysize),datasize(1)),double(shapeid{2}):-1:-double(shapeid{3}),arraysize(1),arraysize(2));
                elseif(strcmpi(shapeid{1},'toeplitz'))
                        arraydata=reshape(arraydata,flipud(datasize(:))');
                        ndata=toeplitz(arraydata(1:arraysize(1),2),arraydata(1:arraysize(2),1));
                end
                if(opt.fullarrayshape && issparse(ndata))
                        ndata=cast(full(ndata),data(j).(N_('_ArrayType_')));
                end
        elseif(isfield(data,N_('_ArraySize_')))
            if(isstruct(data(j).(N_('_ArraySize_'))))
                data(j).(N_('_ArraySize_'))=jdatadecode(data(j).(N_('_ArraySize_')),opt);
            end
            if(iscpx)
                ndata=complex(ndata(1,:),ndata(2,:));
            end
            if(format>1.9 && iscol==0)
                data(j).(N_('_ArraySize_'))=data(j).(N_('_ArraySize_'))(end:-1:1);
            end
            dims=data(j).(N_('_ArraySize_'))(:)';
            if(iscell(dims))
                dims=cell2mat(dims);
            end
	    if(length(dims)==1)
		dims=[1 dims];
	    end
            ndata=reshape(ndata(:),dims(:)');
            if(format>1.9 && iscol==0)
                ndata=permute(ndata,ndims(ndata):-1:1);
            end
        end
        newdata{j}=ndata;
      end
      if(len==1)
          newdata=newdata{1};
      end
    end

    %% handle table data
    if(isfield(data,N_('_TableRecords_')))
        newdata=cell(len,1);
        for j=1:len
            ndata=data(j).(N_('_TableRecords_'));
            if(iscell(ndata))
                if(iscell(ndata{1}))
                    rownum=length(ndata);
                    colnum=length(ndata{1});
                    nd=cell(rownum, colnum);
                    for i1=1:rownum
                        for i2=1:colnum
                            nd{i1,i2}=ndata{i1}{i2};
                        end
                    end
                    newdata{j}=cell2table(nd);
                else
                    newdata{j}=cell2table(ndata);
                end
            else
                newdata{j}=array2table(ndata);
            end
            if(isfield(data(j),N_('_TableRows_'))&& ~isempty(data(j).(N_('_TableRows_'))))
                newdata{j}.Properties.RowNames=data(j).(N_('_TableRows_'))(:);
            end
            if(isfield(data(j),N_('_TableCols_')) && ~isempty(data(j).(N_('_TableCols_'))))
                newdata{j}.Properties.VariableNames=data(j).(N_('_TableCols_'));
            end
        end
        if(len==1)
            newdata=newdata{1};
        end
    end

    %% handle map data
    if(isfield(data,N_('_MapData_')))
        newdata=cell(len,1);
        for j=1:len
            key=cell(1,length(data(j).(N_('_MapData_'))));
            val=cell(size(key));
            for k=1:length(data(j).(N_('_MapData_')))
                key{k}=data(j).(N_('_MapData_')){k}{1};
                val{k}=jdatadecode(data(j).(N_('_MapData_')){k}{2},opt);
            end
            ndata=containers.Map(key,val);
            newdata{j}=ndata;
        end
        if(len==1)
            newdata=newdata{1};
        end
    end

    %% handle graph data
    if(isfield(data,N_('_GraphNodes_')) && exist('graph','file') && exist('digraph','file'))
        newdata=cell(len,1);
        isdirected=1;
        for j=1:len
            nodedata=data(j).(N_('_GraphNodes_'));
            if(isstruct(nodedata))
                nodetable=struct2table(nodedata);
            elseif(isa(nodedata,'containers.Map'))
                nodetable=[keys(nodedata);values(nodedata)];
                if(strcmp(nodedata.KeyType,'char'))
                    nodetable=table(nodetable(1,:)',nodetable(2,:)','VariableNames',{'Name','Data'});
                else
                    nodetable=table(nodetable(2,:)','VariableNames',{'Data'});
                end
            else
                nodetable=table;
            end

            if(isfield(data,N_('_GraphEdges_')))
                edgedata=data(j).(N_('_GraphEdges_'));
            elseif(isfield(data,N_('_GraphEdges0_')))
                edgedata=data(j).(N_('_GraphEdges0_'));
                isdirected=0;
            elseif(isfield(data,N_('_GraphMatrix_')))
                edgedata=jdatadecode(data(j).(N_('_GraphMatrix_')),varargin{:});
            end

            if(exist('edgedata','var'))
                if(iscell(edgedata))
                    endnodes=edgedata(:,1:2);
                    endnodes=reshape([endnodes{:}],size(edgedata,1),2);
                    weight=cell2mat(edgedata(:,3:end));
                    edgetable=table(endnodes,[weight.Weight]','VariableNames',{'EndNodes','Weight'});

                    if(isdirected)
                        newdata{j}=digraph(edgetable,nodetable);
                    else
                        newdata{j}=graph(edgetable,nodetable);
                    end
                elseif(ndims(edgedata)==2 && isstruct(nodetable))
                    newdata{j}=digraph(edgedata,fieldnames(nodetable));
                end
            end
        end
        if(len==1)
            newdata=newdata{1};
        end
    end

    %% handle bytestream and arbitrary matlab objects
    if(isfield(data,N_('_ByteStream_')) && isfield(data,N_('_DataInfo_')))
        newdata=cell(len,1);
        for j=1:len
            if(isfield(data(j).(N_('_DataInfo_')),'MATLABObjectClass'))
                if(needbase64)
                    newdata{j}=getArrayFromByteStream(base64decode(data(j).(N_('_ByteStream_'))));
                else
                    newdata{j}=getArrayFromByteStream(data(j).(N_('_ByteStream_')));
                end
            end
        end
        if(len==1)
            newdata=newdata{1};
        end
    end

    %% handle data link
    if(opt.maxlinklevel>0 && isfield(data,N_('_DataLink_')))
        if(ischar(data.(N_('_DataLink_'))))
            datalink=data.(N_('_DataLink_'));
            if(regexp(datalink,'\:\$'))
                ref=regexp(datalink, '^(?<proto>[a-zA-Z]+://)*(?<path>.+)(?<delim>\:)()*(?<jsonpath>(?<=:)\$\d*\.*.*)*', 'names');
            else
                ref=regexp(datalink, '^(?<proto>[a-zA-Z]+://)*(?<path>.+)(?<delim>\:)*(?<jsonpath>(?<=:)\$\d*\..*)*', 'names');
            end
            if(~isempty(ref.path))
                uripath=[ref.proto ref.path];
                [fpath, fname, fext]=fileparts(uripath);
                opt.maxlinklevel=opt.maxlinklevel-1;
                switch(lower(fext))
                    case {'.json','.jnii','.jdt','.jdat','.jmsh','.jnirs'}
                        newdata=loadjson(uripath, opt);
                    case {'.bjd' ,'.bnii','.jdb','.jbat','.bmsh','.bnirs', '.jamm'}
                        newdata=loadbj(uripath, opt);
                    case {'.ubj'}
                        newdata=loadubjson(uripath, opt);
                    case {'.msgpack'}
                        newdata=loadmsgpack(uripath, opt);
                    case {'.h5','.hdf5','.snirf'}  % this requires EasyH5 toolbox
                        newdata=loadh5(uripath, opt);
                    otherwise
                        % _DataLink_ url does not specify type, assuming JSON format
                        if(regexpi(datalink,'^\s*(http|https|ftp|file)://'))
                            newdata=loadjson(uripath, opt);
                        else
                            warning('_DataLink_ url is not supported')
                        end
                end
                if(~isempty(ref.jsonpath))
                    newdata=getfromjsonpath(newdata,ref.jsonpath);
                end
            end
        end
    end

    %% subfunctions 
    function escaped=N_(str)
      escaped=[prefix str];
    end
end
