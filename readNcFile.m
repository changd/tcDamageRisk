function [long, lat, t, wind, n_tracks] = readNcFile(path)

    file_info = ncinfo(path); % Gets information of the NETCDF file
    VAR = file_info.Variables; % For instance, here are the variables
    S = size(VAR); % Number of variables in your file
    % For example let's try to get the name & size of the variables
    NAME=cell(S); SIZ=zeros(S(2),3);

    for i = 1:S(2)

        NAME(i) = cellstr(VAR(i).Name); % Allocate name of variable

        s = VAR(i).Size; ss = size(s); % s is variable size, ss is number of dim
        SIZ(i,1) = s(1); % Allocate first dim of var
        if ss(2)>=3, SIZ(i,2:3) = s(2:3); % and 2nd & 3rd if 3D
        elseif ss(2)>=2, SIZ(i,2:3) = [s(2) 0]; % 2nd if 2D
        else SIZ(i,2:3) = [0 0]; % Set 2nd & 3rd dim to zero if var is 1D
        end

    end

    dat = SAM_loadVars(NAME(1:end), path);
    long = dat.lon; % u
    lat = dat.lat; % v
    t = dat.t; % w
    wind = dat.wind;
    n_tracks = dat.n_tracks;