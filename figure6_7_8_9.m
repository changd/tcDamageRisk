clearvars; close all; fclose('all');
curdir = pwd; % Current directory
dataPath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/fhlo/Hermine';
% sourcePath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/latex/figures/';

path = [dataPath,'/ecmwf_201609010000_AL092016_lib0', num2str(0), '_DEREK.nc']; % Path of file in current directory
[long, lat, t, wind, n_tracks] = readNcFile(path);

%% Define dimensions // initialize matrices
dimsLong = 200; dimsLat = 100; dimsT = 121; nMembers = 1000; 
nFiles = 10; nMembersPerFile = 100;

% Ensemble forecast is split across 10 files, with 100 ensemble members per
% file. Wind data is indexed by latitude, longitude, time, and ensemble
% member.

% Calculate mean wind field for each file
windCollect = zeros(nMembersPerFile, dimsT, dimsLong, dimsLat);
meanVel = zeros(nFiles, dimsT, dimsLong, dimsLat);
meanVelFunc = zeros(nFiles, dimsT, dimsLong, dimsLat);
meanPoisson = zeros(nFiles, dimsT, dimsLong, dimsLat);

% Calculate CDF (cumulative intensity) for each member
cdf = zeros(nMembers, dimsLong, dimsLat);

%% Get results
for j = 0:9
    % Read ensemble forecast data
    path = [dataPath,'/ecmwf_201609010000_AL092016_lib0', num2str(j), '_DEREK.nc']; % Path of file in current directory
    [long, lat, t, wind, n_tracks] = readNcFile(path);

    % Collect wind data into a large matrix
    for k = 1:100
        windCollect(k,1:dimsT,1:dimsLong,1:dimsLat) = wind(:,k,:,:);
    end
    
    % Get statistics for the set of ensemble members
    tic
    stormStatistics = GetLineFailureModelEnsemble(windCollect);
    toc
    
    meanVel(j+1,:,:,:) = squeeze(mean(windCollect));
    meanVelFunc(j+1,:,:,:) = stormStatistics.windMean;
    meanPoisson(j+1,:,:,:) = stormStatistics.lambdaMean;
    
    % Put cumulative intensities in one large matrix
    idxCdf = j*100+1;
    cdf(idxCdf:idxCdf+99,:,:) = stormStatistics.cLambda;
end

%% Poisson example for one ensemble member of Hermine
ensembleIdx = 1;

% Calculate Poisson intensity for single ensemble member
alpha=4175.6; Vcrit = 20.6; constTerm = 0.62*0.49/365/24;
velFuncPlot = zeros(dimsT, dimsLong, dimsLat);
for i = 1:dimsT
    for j = 1:dimsLong
        for k = 1:dimsLat
            velFuncPlot(i,j,k) = max([wind(i,1,j,k) Vcrit])/Vcrit;
        end
    end
end
poissonPlot = constTerm*(1-alpha)*ones(size(velFuncPlot)) + ...
    constTerm*alpha*(velFuncPlot).^2;

figure
% Plot velocities
for i = 1:6
    timeStep = 27+2*i;
    
    subplot(2,6,i)
    [C,h] = contourf(long, lat, squeeze(wind(timeStep,ensembleIdx,:,:))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 45])
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    
    if timeStep-24 < 10
        title(['2016/09/02, 0' num2str(timeStep-24) 'z'], 'Interpreter', 'latex')        
    else
        title(['2016/09/02, ' num2str(timeStep-24) 'z'], 'Interpreter', 'latex')
    end
end

% Plot Poisson intensities
for i = 1:6
    timeStep = 27+2*i;
    
    subplot(2,6,i+6)
    [C,h] = contourf(long, lat, 3.4*squeeze(poissonPlot(timeStep,:,:))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 2])
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1800 600])

%% Plot ensemble-averaged mean velocities 
meanVelPlot = squeeze(mean(meanVel)); 

figure
for i = 1:6
    timeStep = 27+2*i;

    subplot(1,6,i)
    [C,h] = contourf(long, lat, squeeze(meanVelPlot(timeStep,:,:))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 45])
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    
    if timeStep-24 < 10
        title(['2016/09/02, 0' num2str(timeStep-24) 'z'], 'Interpreter', 'latex')        
    else
        title(['2016/09/02, ' num2str(timeStep-24) 'z'], 'Interpreter', 'latex')
    end
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1800 300])

%% Calculate cumulative intensity as a function of mean wind field
meanVelFuncPlot = zeros(size(meanVelPlot));
for i = 1:dimsT
    for j = 1:dimsLong
        for k = 1:dimsLat
            meanVelFuncPlot(i,j,k) = max([meanVelPlot(i,j,k) Vcrit])/Vcrit;
        end
    end
end
meanPoissonPlot = constTerm*(1-alpha)*ones(size(meanVelPlot)) + ...
    constTerm*alpha*(meanVelFuncPlot).^2;
cdfMeanVel = squeeze(sum(meanPoissonPlot));

lengthLine = 7.08;

figure
subplot(1,2,1)
[C,h] = contourf(long, lat, lengthLine*cdfMeanVel');
caxis([0 4.5]); c = colorbar; c.TickLabelInterpreter = 'latex';
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) max(long)]); ylim([min(lat) max(lat)]);
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)

subplot(1,2,2)
[C,h] = contourf(long, lat, lengthLine*squeeze(mean(cdf))');
caxis([0 4.5]); c = colorbar; c.TickLabelInterpreter = 'latex';
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) max(long)]); ylim([min(lat) max(lat)]);
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)

set(gcf,'Renderer', 'painters', 'Position', [10 10 800 260])

%% Plot cumulative intensities in subregion 
% Indicates locations of where we conduct further analysis
lengthLine = 7.08;
idxLat = find(lat <= 29.905 & lat >= 29.895);
cdfPlot = lengthLine*cdf(:,:,idxLat);

% Calculate number of Lambdas (cumulative intensities) in each bin
binBounds = 0:0.1:10; bins = zeros(length(binBounds)-1, length(long));
for i = 1:length(binBounds)-1
    for j = 1:length(long)
        idx = find(cdfPlot(:,j) >= binBounds(i) & cdfPlot(:,j) < binBounds(i+1));
        bins(i,j) = length(idx);
    end
end

%%

% Plot CDF
figure
[C,h] = contourf(long, lat, lengthLine*squeeze(mean(cdf))');

c = colorbar; c.TickLabelInterpreter = 'latex';
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) -83.8]); ylim([min(lat) 30.1]);

set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)
hold on

% Indicate locations (idxPlot) where we will plot histograms of Lambda
plot([min(long); max(long)], lat(idxLat)*ones(2,1), '-k', 'LineWidth', 2)
hold on

idxPlot = 1:40:40*3+1;
plot(long(idxPlot), lat(idxLat)*ones(1,length(idxPlot)), '.k', 'MarkerSize', 24)
xticks(-85.2:0.1:-83.8)

set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 200])

%% Plot histogram of Lambdas at select locations (Figure 20)
figure
for i = 1:4
    idxPlot = 1+(i-1)*40;
    subplot(1,4,i)
    histogram(cdfPlot(:,idxPlot), 0:0.5:10)
    xlabel('Failure Rate (hr$^{-1}$)','Interpreter','latex'); 
    xlim([min(binBounds) max(binBounds)]); ylim([0 600]);
    xticks([0:2:10]); yticks([0:150:600])
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 300])

%% Calculate failure probabilities under different models
noFailsMax = 40; cdfMatLine = lengthLine*cdf;

% FD-A
cdfMeanPoisson = mean(cdf);
probFailMean = zeros(noFailsMax, length(long), length(lat));
cdfMeanLine = lengthLine*cdfMeanPoisson;
for i = 1:noFailsMax+1
    probFailMini = cdfMeanLine.^(i-1) .* exp(-cdfMeanLine) / factorial(i-1);
    probFailMean(i,:,:) = probFailMini;
end

for i = 1:length(long)
    for j = 1:length(lat)
        probFailMean(:,i,j) = probFailMean(:,i,j)/sum(probFailMean(:,i,j));
    end
end

% FD-B
probFail = zeros(noFailsMax, length(long), length(lat));
for i = 1:noFailsMax+1
    probFailMini = cdfMatLine.^(i-1) .* exp(-cdfMatLine) / factorial(i-1);
    probFailMiniAvg = mean(probFailMini);
    probFail(i,:,:) = probFailMiniAvg;
end

for i = 1:length(long)
    for j = 1:length(lat)
        probFail(:,i,j) = probFail(:,i,j)/sum(probFail(:,i,j));
    end
end

%% Plot failure probabilities (Figure 20)
probFailMeanAnalyze = zeros(9, 4); probFailAnalyze = zeros(9, 4); 
probExceedMean = zeros(4, 1); probExceed = zeros(4, 1);
noExceedMean = zeros(4, 1); noExceed = zeros(4, 1);

figure
for i = 1:4
    idxPlot = 1+(i-1)*40;

    subplot(1,4,i)
    plot(0:noFailsMax, probFailMean(:, idxPlot, idxLat), 'LineWidth', 2)
    hold on
    plot(0:noFailsMax, probFail(:, idxPlot, idxLat), 'LineWidth', 2)
    xlabel('\# of Failures (hr$^{-1}$)','Interpreter','latex'); ylabel('Probability','Interpreter','latex')
    xlim([0 18]); ylim([0 0.6])
    xticks(0:3:18)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    legend({'FD-A', 'FD-B'}, 'Interpreter', 'latex')
    
    probFailMeanAnalyze(:, i) = probFailMean(1:9, idxPlot, idxLat);
    probFailAnalyze(:, i) = probFail(1:9, idxPlot, idxLat);
    
    probExceedMean(i) = sum(probFailMean(10:end, idxPlot, idxLat));
    probExceed(i) = sum(probFail(10:end, idxPlot, idxLat));
    
    probTotal = 0; iter = 0;
    while probTotal < 0.95
        iter = iter + 1;
        probTotal = probTotal + probFailMean(iter, idxPlot, idxLat);
    end
    noExceedMean(i) = iter - 2;
    
    probTotal = 0; iter = 0;
    while probTotal < 0.95
        iter = iter + 1;
        probTotal = probTotal + probFail(iter, idxPlot, idxLat);
    end
    noExceed(i) = iter - 1;
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 300])
% saveas(gcf, [sourcePath, 'hermineNoFailures.jpg'])