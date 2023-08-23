function recursiveDownload(weblocation, localFolder)
    
    % Read the HTML content of the URL
    htmlContent = webread(weblocation);    
    pattern = '<a href="([^"]+)">';
    matches = regexp(htmlContent, pattern, 'tokens');
    
    % Iterate over the matches
    for i = 2:numel(matches) % Ignore i=1, which is the parent directory link: '../'
        item = matches{i}{1};
        
        if endsWith(item, '/') % Folder
            % Create the necessary directories if they do not exist
            subfolder = fullfile(localFolder, item);
            if ~isfolder(subfolder)
                mkdir(subfolder);
            end
            
            % Recursively download the subfolder
            subWebLocation = strcat(weblocation, '/', item);
            recursiveDownload(subWebLocation, subfolder);
        else % File
            % Create the necessary directories if they do not exist
            
            if ~isfolder(localFolder)
                mkdir(localFolder);
            end

            % Download the file
            fileUrl = strcat(weblocation, '/', item);
            localFilePath = fullfile(localFolder, item);
            websave(localFilePath, fileUrl);
        end
    end
end
