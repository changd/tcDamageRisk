function stormStatistics = GetLineFailureModelEnsemble(wind) 

    %% Get storm statistics for an ensemble forecast
    % Input: wind -- 4D wind field data (indexed by ensemble, time, lat/long)
    % Output: stormStatistics -- wind velocities, Poisson intensities,
    % cumulative intensities

    % Obtain Poisson intensities, cumulative intensities
        alpha=4175.6; Vcrit = 20.6; constTerm = 0.62*0.49/365/24;
    parameters = [constTerm alpha Vcrit];
    
    % Calculate velocity Statistics
    [wind, windMean] = getVelocityStatistics(wind, Vcrit);
    
    % Calculate Poisson intensities
    [lambda, lambdaMean, lambdaMeanWind] ... 
        = getPoissonIntensityStatistics(wind, windMean, parameters);
    
    % Calculate cumulative intensities
    [cLambda, cLambdaMean, cLambdaMeanWind] ...
        = getCumulativeIntensityStatistics(lambda, lambdaMean, lambdaMeanWind);
    
    stormStatistics = struct([]);
    stormStatistics(1).wind = wind; 
    stormStatistics(1).windMean = windMean;
    stormStatistics(1).lambda = lambda;
    stormStatistics(1).lambdaMean = lambdaMean;
    stormStatistics(1).lambdaMeanWind = lambdaMeanWind;
    stormStatistics(1).cLambda = cLambda;
    stormStatistics(1).cLambdaMean = cLambdaMean;
    stormStatistics(1).cLambdaMeanWind = cLambdaMeanWind;
end

function [fV, windMean] = getVelocityStatistics(wind, Vcrit)
    fV = max(Vcrit*ones(size(wind)), wind);
    windMean = squeeze(mean(wind));    
end

function [lambda, lambdaMean, lambdaMeanWind] ...
    = getPoissonIntensityStatistics(wind, saaWindMean, parameters)

    constTerm = parameters(1); alpha = parameters(2); Vcrit = parameters(3);
    Poisson = @(V, alpha, Vcrit) ones(size(V)) + alpha*((V/Vcrit).^2 - ones(size(V)));    

    % Calculate Poisson intensities
    lambda = constTerm * Poisson(wind, alpha, Vcrit);
    
    % Calculate Poisson intensity averaged across ensemble members
    lambdaMean = squeeze(mean(lambda));
    
    % Calculate Poisson intensity as function of mean velocity
    lambdaMeanWind = constTerm * Poisson(saaWindMean, alpha, Vcrit);
end

function [cLambda, cLambdaMean, cLambdaMeanWind] ...
    = getCumulativeIntensityStatistics(lambda, lambdaMean, lambdaMeanWind)

    % Calculate cumulative Poisson intensities
    cLambda = squeeze(sum(lambda, 2));
    
    % Calculuate cumulative Poisson intensities averaged across ensemble
    % members
    cLambdaMean = squeeze(sum(lambdaMean)); 
    
    % Calculate cumulative Poisson intensity as function of mean velocity
    cLambdaMeanWind = squeeze(sum(lambdaMeanWind));
end