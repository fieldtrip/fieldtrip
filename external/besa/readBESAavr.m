function avr = ReadBESAavr(filename)

% Reads BESA *.avr files

fp = fopen(filename,'r');
if (fp)
    % Read header of .avr file
    avr.Npts = fscanf(fp,'Npts= %i');
    avr.TSB = fscanf(fp,' TSB= %f');    
    avr.DI = fscanf(fp,'DI= %f');
    avr.Time = [avr.TSB:avr.DI:avr.TSB+avr.DI*(avr.Npts-1)];
    temp = fscanf(fp,' SB= %f SC= %f');
    try 
        avr.Nchan = fscanf(fp,'Nchan= %i\n'); 
    catch
    end

    if avr.Nchan >0
        avr.ChannelLabels = fgetl(fp);
    end

    for i=1:1000
        try avr.Data(i,:)=fscanf(fp,'%f',[avr.Npts,1]);
        catch
            break
        end
    end
    avr.Nchan=size(avr.Data,1);
end
fclose(fp);
