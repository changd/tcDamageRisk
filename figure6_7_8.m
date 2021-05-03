clearvars; close all; fclose('all');
curdir = pwd; % Current directory
dataPath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/fhlo/Hermine';
sourcePath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/latex/figures/';

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

% Plot CDF
figure
[C,h] = contourf(long, lat, lengthLine*squeeze(mean(cdf))');
caxis([1 4.5]); c = colorbar; c.TickLabelInterpreter = 'latex';
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) -83.8]); ylim([min(lat) 30.1]);
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)
hold on

% Indicate locations (idxPlot) where we will plot histograms of Lambda
plot([min(long); max(long)], lat(idxLat)*ones(2,1), '-k', 'LineWidth', 2)
hold on
idxPlot = 1:25:25*5+1;
xticks(-85.2:0.1:-83.8)
plot(long(idxPlot), lat(idxLat)*ones(1,6), '.k', 'MarkerSize', 24)
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 200])

%% Plot histogram of Lambdas at select locations 
figure
for i = 1:6
    idxPlot = 1+(i-1)*25;
    subplot(1,6,i)
    histogram(cdfPlot(:,idxPlot), 0:0.5:10)
    xlabel('Failure Rate','Interpreter','latex'); 
    xlim([min(binBounds) max(binBounds)]); ylim([0 600]);
    xticks(0:2:8); yticks(0:150:600)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 200])

%% Calculate failure probabilities under different models
noFailsMax = 40; cdfMatLine = 3.4*cdf;

% Model A
cdfMeanPoisson = mean(cdf);
probFailMean = zeros(noFailsMax, length(long), length(lat));
cdfMeanLine = 3.4*cdfMeanPoisson;
for i = 1:noFailsMax+1
    probFailMini = cdfMeanLine.^(i-1) .* exp(-cdfMeanLine) / factorial(i-1);
    probFailMean(i,:,:) = probFailMini;
end

for i = 1:length(long)
    for j = 1:length(lat)
        probFailMean(:,i,j) = probFailMean(:,i,j)/sum(probFailMean(:,i,j));
    end
end

% Model B
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
figure
for i = 1:6
    idxPlot = 1+(i-1)*25;

    subplot(1,6,i)
    plot(0:noFailsMax, probFailMean(:, idxPlot, idxLat), 'LineWidth', 2)
    hold on
    plot(0:noFailsMax, probFail(:, idxPlot, idxLat), 'LineWidth', 2)
    xlabel('No. Failures (hr$^{-1}$km$^{-1}$','Interpreter','latex'); ylabel('Probability','Interpreter','latex')
    xlim([0 8]); ylim([0 0.6])
    xticks(0:2:8)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    legend({'Model A', 'Model B'}, 'Interpreter', 'latex')
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1300 200])
saveas(gcf, [sourcePath, 'herminePlotNoFailures.jpg'])

%%
idxLat = find(lat <= 30.605 & lat >= 30.595);
cdfPlot = 3.4*cdf(:,:,idxLat);

binBounds = 0:0.1:10; bins = zeros(length(binBounds)-1, length(long));
for i = 1:length(binBounds)-1
    for j = 1:length(long)
        idx = find(cdfPlot(:,j) >= binBounds(i) & cdfPlot(:,j) < binBounds(i+1));
        bins(i,j) = length(idx);
    end
end

%%
figure

imagesc(long, binBounds, bins)
c = colorbar;
c.TickLabelInterpreter = 'latex';
title('Empirical Probability of $\Lambda$', 'Interpreter', 'latex')
xlabel('Longitude','Interpreter','latex'); ylabel('$\Lambda$','Interpreter','latex')
% xlim([min(distPlot(100:idxStop)) max(distPlot(100:idxStop))]); 
ylim([min(binBounds) 0.5]);
set(gca, 'YDir', 'normal')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)
% saveas(gcf, [sourcePath, 'cdfImageSC.jpg'])

%%
figure

for i = 1:8
    subplot(2,4,i)
    histogram(cdfPlot(:,i*25))
    % title('$r = 40$ km', 'Interpreter', 'latex')
    xlabel('$\Lambda$','Interpreter','latex'); %ylabel('$\Lambda$','Interpreter','latex')
    xlim([min(binBounds) max(binBounds)]); ylim([0 800]);
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end

%%
figure

subplot(1,2,1)
imagesc(long, 0:noFailsMax, probFailMean(:,:,idxLat))
c = colorbar;
c.TickLabelInterpreter = 'latex';
% caxis([0 0.8])
title('Model A', 'Interpreter', 'latex')
xlabel('Longitude','Interpreter','latex'); ylabel('Number of failures','Interpreter','latex')
% xlim([min(distPlot(100:idxStop)) max(distPlot(100:idxStop))]); 
ylim([0 8]);
set(gca, 'YDir', 'normal')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)

subplot(1,2,2)
imagesc(long, 0:noFailsMax, probFail(:,:,idxLat))
c = colorbar;
c.TickLabelInterpreter = 'latex';
% caxis([0 0.8])
title('Model B', 'Interpreter', 'latex')
xlabel('Longitude','Interpreter','latex'); ylabel('Number of failures','Interpreter','latex')
% xlim([min(distPlot(100:idxStop)) max(distPlot(100:idxStop))]); 
ylim([0 8]);
set(gca, 'YDir', 'normal')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)
set(gcf,'Renderer', 'painters', 'Position', [10 10 800 400])
% saveas(gcf, [sourcePath, 'probFailImageSC.jpg'])

%% Calculate crit zone size
idx = find(cdfMeanVel > 0.0042);
FR1 = lengthLine*cdfMeanVel;
FR1_mean = mean(mean(cdfMeanVel(idx)));

idx = find(mean(cdf) > 0.0042);
FR2 = lengthLine*squeeze(mean(cdf));
FR2_mean = mean(mean(FR2(idx)))/lengthLine;

%% Calculate velocity exceedances
idxExceed = zeros(dimsT, 1);
meanExceedance = zeros(dimsT, 1);
for i = 1:dimsT
    idxExceed(i) = length(find(meanVelPlot(i,:,:) > 20.6));
    idxs = find(meanVelPlot(i,:,:) > 20.6);
    if length(idxExceed(i)) > 0
        velExceed = meanVelPlot(i, idxs) - 20.6;
        meanExceedance(i) = mean(mean(velExceed));
    end
end

