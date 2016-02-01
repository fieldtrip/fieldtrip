function simulink2svg
% Show the hidden handles to get access to the simulink figures
set(0,'ShowHiddenHandles','on');
% Iterate over all children of the root
child = get(0,'children');                                
for i = 1:length(child)                                  
    % Pick all scope figures
    if strcmp(get(child(i),'Tag'),'SIMULINK_SIMSCOPE_FIGURE')
        % In order to get the correct background turn inverted background
        % off. If you like a white background you should invert all labels
        % and lines.
        set(child(i),'InvertHardcopy','off');
        % Use plot2svg
        plot2svg(['Scope_' num2str(i) '.svg'], child(i))          
        set(0,'ShowHiddenHandles','off')                         
    end                                                      
end                                       