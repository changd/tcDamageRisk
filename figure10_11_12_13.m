clearvars; close all; fclose('all');
curdir = pwd; % Current directory
dataPath = '/Users/changd/Dropbox (MIT)/UROP_2020/urop2020_Neel/fhlo_michael/';
sourcePath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/latex/ress/figures/';

%% Define dimensions // initialize matrices
%can be used to downsize dimensions
start_long = -86.5;
end_long = -83.4;
start_lat = 29.6;
end_lat = 32.2;

%set from readNcFile
long_nc = -86.5:0.01:-83.4;
lat_nc = 29.6:0.01:32.2;


%our downsized arrays
idx_long = find(long_nc <= (end_long+ 0.001) & long_nc >= (start_long-0.001));
idx_lat = find(lat_nc <= (end_lat+ 0.001) & lat_nc >= (start_lat-0.001));  %stupid floating point issues
long = start_long:0.01:end_long;
lat = start_lat:0.01:end_lat;
long_offset = idx_long(1);
lat_offset = idx_lat(1);


dimsLong = size(long, 2); 
dimsLat = size(lat, 2);
dimsT = 121; 
nFiles = 10; 
nMembersPerFile = 50;
nMembers = nFiles * nMembersPerFile;

%% Poisson example for one ensemble member of Michael (Figure 17)
load('/Users/changd/Dropbox (MIT)/UROP_2020/urop2020_sharedMaterialsDerek/michaelData_v2.mat');

ensembleIdx = 1;

% Calculate Poisson intensity for single ensemble member
alpha=4175.6; Vcrit = 20.6; constTerm = 0.62*0.49/365/24;
velFuncPlot = zeros(dimsT, dimsLong, dimsLat);
for i = 1:dimsT
    i
    for j = 1:dimsLong
        for k = 1:dimsLat
            velFuncPlot(i,j,k) = max([wind(i,1,j+long_offset - 1,k+lat_offset - 1) Vcrit])/Vcrit;
        end
    end
end
poissonPlot = constTerm*(1-alpha)*ones(size(velFuncPlot)) + ...
    constTerm*alpha*(velFuncPlot).^2;

%%
lengthLine = 7.08;

figure
% Plot velocities
for i = 1:6
    timeStep = 29+2*i;
    
    subplot(2,6,i)
    [C,h] = contourf(long, lat, squeeze(wind(timeStep,ensembleIdx,idx_long,idx_lat))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 54]);
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    if timeStep - 24 < 10
        title(['2020/10/11, 0' num2str(timeStep - 24), 'z'], 'Interpreter', 'latex')
    elseif timeStep - 24 >= 10
        title(['2020/10/11, ' num2str(timeStep-24), 'z'], 'Interpreter', 'latex')
    end
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end

% Plot Poisson intensities
for i = 1:6
    timeStep = 29+2*i;
    
    subplot(2,6,i+6)
    [C,h] = contourf(long, lat, lengthLine * squeeze(poissonPlot(timeStep,:,:))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 5])
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1800 600])

%% Plot cumulative intensity for single ensemble member
figure
[C,h] = contourf(long, lat, 3.4*squeeze(sum(poissonPlot))');
c = colorbar; c.TickLabelInterpreter = 'latex';
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)

%%
meanVelPlot = squeeze(mean(meanVel)); 

figure
for i = 1:6
    timeStep = 29+2*i;

    subplot(1,6,i)
    [C,h] = contourf(long, lat, squeeze(meanVelPlot(timeStep,:,:))');
    c = colorbar; c.TickLabelInterpreter = 'latex'; caxis([0 54])
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    if timeStep - 24 < 10
        title(['2020/10/11, 0' num2str(timeStep - 24), 'z'], 'Interpreter', 'latex')
    elseif timeStep - 24 >= 10
        title(['2020/10/11, ' num2str(timeStep-24), 'z'], 'Interpreter', 'latex')
    end    
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
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

%% Plot comparison of cumulative intensities 
lengthLine = 7.08;

figure
subplot(1,2,1)
[~,h] = contourf(long, lat, lengthLine*cdfMeanVel');
caxis([0 14]); 
c = colorbar; c.TickLabelInterpreter = 'latex';
% title('Model 2', 'Interpreter', 'latex')
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) max(long)]); ylim([min(lat) max(lat)]);
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 18)

subplot(1,2,2)
[C,h] = contourf(long, lat, lengthLine*squeeze(mean(cdf))');
hold on
caxis([0 14]); 
c = colorbar; c.TickLabelInterpreter = 'latex';
% img = imread([dataPath 'county_map.jpg']);
% county_map = imshow(img);
% set(county_map, 'AlphaData', 0.75);
% image('CData', img);
% title('Model 3', 'Interpreter', 'latex')
xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([min(long) max(long)]); ylim([min(lat) max(lat)]);
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 18)

set(gcf,'Renderer', 'painters', 'Position', [10 10 900 300])

%% Plot cumulative intensities in subregion
% (indicates locations of where we will conduct further anaysis)
idxLat = find(lat <= 29.755 & lat >= 29.745);
cdfPlot = lengthLine*cdf(:,:,idxLat);

% Calculate number of Lambdas (cumulative intensities) in each bin
binBounds = 0:0.1:14; bins = zeros(length(binBounds)-1, length(long));
for i = 1:length(binBounds)-1
    for j = 1:length(long)
        idx = find(cdfPlot(:,j) >= binBounds(i) & cdfPlot(:,j) < binBounds(i+1));
        bins(i,j) = length(idx);
    end
end

% Plot CDF
figure
[C,h] = contourf(long, lat, lengthLine*squeeze(mean(cdf))');
%caxis([1 2.2]); 

c = colorbar; c.TickLabelInterpreter = 'latex';

xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
xlim([-86.5 -85.6]); ylim([29.6 31]);

set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 16)
hold on

% Indicate locations (idxPlot) where we will plot histograms of Lambda
plot([min(long); max(long)], lat(idxLat)*ones(2,1), '-k', 'LineWidth', 2)
hold on

idxPlot = 5:27:5+27*3;
% xticks(start_long:0.1:end_long)
plot(long(idxPlot), lat(idxLat)*ones(1,length(idxPlot)), '.k', 'MarkerSize', 24)

set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 200])

%% Plot histogram of Lambdas at select locations 
figure
for i = 1:4
    idxPlot = 5+(i-1)*27;
    
    subplot(1,4,i)
    histogram(cdfPlot(:,idxPlot), 0:0.5:14)
    
    xlabel('Failure Rate (hr$^{-1}$)','Interpreter','latex'); 
    xlim([min(binBounds) max(binBounds)]); ylim([0 200]);
    xticks(0:2:14); yticks(0:50:200)
    
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
figure
for i = 1:4
    idxPlot = 1+(i-1)*27;

    subplot(1,4,i)
    plot(0:noFailsMax, probFailMean(:, idxPlot, idxLat), 'LineWidth', 2)
    hold on
    plot(0:noFailsMax, probFail(:, idxPlot, idxLat), 'LineWidth', 2)
    xlabel('\# of Failures (hr$^{-1}$)','Interpreter','latex'); ylabel('Probability','Interpreter','latex')
%     title(['$29.75^{\circ}N$ ' num2str(round(long(idxPlot), 2)) '$^{\circ}W$'], 'Interpreter', 'latex')
    xlim([0 30]); ylim([0 0.4])
    xticks(0:5:30)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    legend({'FD-A', 'FD-B'}, 'Interpreter', 'latex')
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 1200 300])

%% Calculate crit zone size
idx = find(cdfMeanVel > 0.0042);
length(idx)
FR1 = lengthLine*cdfMeanVel;
FR1_mean = mean(mean(FR1));

idx = find(mean(cdf) > 0.0042);
length(idx)
FR2 = lengthLine*squeeze(mean(cdf));
FR2_mean = mean(mean(FR2));

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
