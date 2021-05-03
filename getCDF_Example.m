function [windField, lambda, CDF] = getCDF_Example(parameters)

    %% Input: set of hurricane/track parameters
    %% Output:
    % windField -- wind velocities at all denoted locations/times
    % lambda -- Poisson intensities at all denoted locations/times
    % CDF -- cumulative Poisson intensities at all denoted locations

    % Get hurricane parameters
    Vtrans = parameters.Vtrans; Vm = parameters.Vm; Rm = parameters.Rm; B = parameters.B;
    
    % Get geographical parameters
    longInit = parameters.longInit; latInit = parameters.latInit; % Storm genesis location
    long = parameters.long; lat = parameters.lat; % lat/long at which velocities are estimated
    latGrid = parameters.latGrid; longGrid = parameters.longGrid; % gridded lat/long locations
    latTrack = parameters.latTrack; longTrack = parameters.longTrack;
    
    % Get storm track locations
    nTimes = length(parameters.latTrack);
    
    %% Get wind fields at each time step
    nLat = length(lat); nLong = length(long);
    
    % Function handle for Holland model-estimated wind velocities
    V_Hol = @(Vm, Rm, B, R) Vm * ((Rm./R).^B .* exp(ones(size(R)) - (Rm./R).^B)).^0.5;

    % Calculate wind field by iterating through all time steps
    windField = zeros(nTimes, nLong, nLat);
    for i = 1:nTimes
        R = pos2distVec(latTrack(i), longTrack(i), latGrid, longGrid, 2);
        windField(i,:,:) = V_Hol(Vm, Rm, B, R);
    end

    %% Get Poisson intensities
    lambda = zeros(nTimes, nLong, nLat);
    
    % Poisson model parameters
    alpha=4175.6; Vcrit = 20.6; lambdaNorm = 0.62*0.49/365/24;
    
    % Poisson intensity function handle
    Poisson = @(V, alpha, Vcrit) ones(size(V)) + alpha*((V/Vcrit).^2 - ones(size(V))); 

    for i = 1:nTimes
        fWindField = max(Vcrit*ones(size(windField(i,:,:))), windField(i,:,:));
        lambda(i,:,:) = lambdaNorm * Poisson(fWindField, alpha, Vcrit);
    end

    %% Get cumulative Poisson intensities
    CDF = squeeze(sum(lambda));
