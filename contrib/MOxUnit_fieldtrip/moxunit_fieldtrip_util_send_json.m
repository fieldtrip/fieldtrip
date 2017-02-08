function result=moxunit_fieldtrip_util_send_json(url, content)
% Send JSON content to URL, or see which methods are available
%
% methods=moxunit_fieldtrip_util_send_json();
% result==moxunit_fieldtrip_util_send_json(url,content);
%
% Inputs:
%   url                 URL to which content is to be sent.
%   content             Struct with JSON content to be sent.
%
% Output:
%   methods             If no inputs are used, then the output is a cell
%                       string with available methods; possible values
%                       are 'webwrite' and 'curl'.
%

    [keys,send_funcs]=get_methods();

    switch nargin

        case 0
            result=keys;

        case 2
            if isempty(send_funcs)
                error('No method available to send output');
            end
            first_send_func=send_funcs{1};
            first_send_func(url, content);

        otherwise
            error('0 or 2 input arguments are required');
    end


function [keys,send_func]=get_methods()
    prop=get_all_potential_methods();


    % see which ones are availabel
    all_keys=fieldnames(prop);
    keep=cellfun(@(x)prop.(x).available(),all_keys);

    keys=all_keys(keep);
    send_func=cellfun(@(x)prop.(x).func,keys,...
                    'UniformOutput',false);


function prop=get_all_potential_methods()
    prop=struct();

    % recent Matlab's webwrite
    prop.webwrite.available=@()~isempty(which('webwrite')) && ...
                                ~isempty(which('weboptions'));
    prop.webwrite.func=@matlabge2015a_webwrite_send;

    % Unix curl (Octave and older Matlab)
    prop.curl.available=@()isunix() && unix('which curl>/dev/null')==0;
    prop.curl.func=@matlablt2015a_webwrite_send;

function matlabge2015a_webwrite_send(url, content)
    options = weboptions('MediaType','application/json');
    webwrite(url, content, options);


function matlablt2015a_webwrite_send(varargin)
    % use webwrite compatibility function (based on curl command line)
    [unused,ft_root_dir]=ft_version();
    compat_dir=fullfile(ft_root_dir,'compat','matlablt2015a');

    % go to compat directory, ensuring that webwrite is used there
    orig_pwd=pwd();
    cleaner=onCleanup(@()cd(orig_pwd));
    cd(compat_dir);

    % run webwrite
    webwrite(varargin{:});

