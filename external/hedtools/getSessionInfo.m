function [cookie, csrftoken] = getSessionInfo(csrf_url)
%% Setup the session for accessing the HED webservices
%  Parameters:
%       csrf_url  = URL for the services
%
%  Returns: 
%        cookie = a string cookie value
%        csrftoken = a string csrf token for the session.
%
    request = matlab.net.http.RequestMessage;
    uri = matlab.net.URI(csrf_url);
    response1 = send(request,uri);
    cookies = response1.getFields('Set-Cookie');
    cookie = cookies.Value;
    data = response1.Body.char;
    csrfIdx = strfind(data,'csrf_token');
    tmp = data(csrfIdx(1)+length('csrf_token')+1:end);
    csrftoken = regexp(tmp,'".*?"','match');
    csrftoken = string(csrftoken{1}(2:end-1));
