clear

%% Parameters
parameters = struct([]);
parameters(1).longInit = -66.45; 
parameters.latInit = 18.25 - 60*0.1; parameters.latFinal = 18.25 + 60*0.1;
parameters.lat = 11:0.01:25; parameters.long = -70:0.01:-63;
[parameters.latGrid, parameters.longGrid] = meshgrid(parameters.lat, parameters.long);

% Hurricane parameters
Vtrans = 0.1; Vm = 37; Rm = [20; 30; 40]; B = 1; 
parameters.Vtrans = Vtrans; 
parameters.Vm = Vm; parameters.B = B;

% Hurricane track lat/long points
parameters.latTrack = parameters.latInit:Vtrans:parameters.latFinal;
nTimes = length(parameters.latTrack);
parameters.longTrack = parameters.longInit*ones(1,nTimes);

%%
failureDistResults = struct([]); 
for i = 1:length(Rm)
    i
    parameters.Rm = Rm(i);
    [failureDistResults(i).windField, failureDistResults(i).lambda, ...
        failureDistResults(i).CDF] = getCDF_Example(parameters);
end

%%
lat = parameters.lat; long = parameters.long;
categories = ["Tropical Storm", "Category 1", "Category 2"];

figure
for i = 1:length(Rm)
    subplot(1,3,i)
    [C,h] = contourf(long, lat, failureDistResults(i).CDF');
    caxis([0 14])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    title(['$R_m$ = ' num2str(Rm(i))], 'Interpreter', 'latex')
    xlabel('Longitude','Interpreter','latex'); ylabel('Latitude','Interpreter','latex')
    xlim([min(long) max(long)]); ylim([min(lat) max(lat)]);
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
end
set(gcf,'Renderer', 'painters', 'Position', [10 10 800 600])

%%
% lambdaNorm = 0.62*0.49/365/24;
% critZoneSize = zeros(3, 1);
% maxFR = zeros(3, 1); 
% totalFailures = zeros(3,1);
% for i = 1:3
%     margin = lambdaNorm*nTimes/1000;
%     idx = find(failureDistResults(i).CDF > lambdaNorm*nTimes + margin);
%     critZoneSize(i) = length(idx);
%     maxFR(i) = max(max(failureDistResults(i).CDF));
%     totalFailures(i) = mean(mean(failureDistResults(i).CDF(idx)));
% end
