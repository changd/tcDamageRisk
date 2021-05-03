  %% groundTruthing_visualization.m
% Plots county-based failure rates and outage percentages (visual plot)

clear

%% Get county information
% Get coordinates of each county's bounding box
S = shaperead('outageData/Florida_Counties.shp', 'usegeocoords', true);

% Get indices of counties that we analyze
needed_counties = readmatrix('outageData/shapefilecountynames', 'Range', 'A:A');

% Get only the counties we need, in alphabetical order
for j = 1:length (needed_counties)
    Counties(j) = S(needed_counties(j));
end

%% Get areas for each county
areaShpFile = readmatrix(['outageData/Guam Data CSV Files/excel/Power_Outages_2018', ...
        convertStringsToChars("1010_1635")], 'Range', 'K:K');
areaAlachua = 2509.7;
convFactor = areaShpFile(1)/areaAlachua;
area = areaShpFile/convFactor;

%% Times for the calculations
calc_times = [34 38 61 64 88];  %these are the labels for the times in the cdf files
outage_times = ["1010_1635", "1010_1950", "1011_1940", "1011_2200", "1012_2305"];  %times in the power outage files
real_times = ["10/10 4:35 PM", "10/10 7:50 PM", "10/11 7:40 PM", "10/11 10:00 PM", "10/12 11:05 PM"];
cnty_cdf_Collect = zeros(37, 5);

for time_step = 1:5
    cnty_cdf = readmatrix(['outageData/County_CDF_Calculations/county_cdf_t=' num2str(calc_times(time_step))])';
    cnty_cdf_Collect(:,time_step) = cnty_cdf;
    
    % Total v. Velocity
    
    pct_outages = 100*readmatrix(['outageData/Guam Data CSV Files/excel/Power_Outages_2018', ...
        convertStringsToChars(outage_times(time_step))], 'Range', 'F:F');
    
    %Color in the counties
    figure
    subplot(2,1,1);
    
    geoshow(Counties)
    for j = 1:length(Counties)
        x = Counties(j).Lon;
        y = Counties(j).Lat;   % Picked the 3rd region
        idx = find(isnan(x)) ;     % find positions of NaNs
        idx = [1 idx length(x)] ;  % append first and last position
        hold on
        for i = 1:length(idx)-1
            pos = idx(i):idx(i+1) ;
            xi = x(pos) ; yi = y(pos) ;
            % Remove NaN's
            xi(isnan(xi)) = [] ;
            yi(isnan(yi)) = [] ;
            patch(xi,yi,pct_outages(j)) ;     % color the area
%             patch(xi,yi,pct_outages_per_area(j)) ;     % color the area
            c = colorbar;
            c.TickLabelInterpreter = 'latex';
            caxis([0 100])
            set(c, 'YTick', 0:20:100)
        end
    end
    xlabel('Longitude', 'Interpreter', 'latex'); ylabel('Latitude', 'Interpreter', 'latex')
    title('Outages per 100 Households', 'Interpreter', 'latex', 'FontSize', 16)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)    
    
    subplot(2,1,2);
    geoshow(Counties)
    for j = 1:length(Counties)
        x = Counties(j).Lon;
        y = Counties(j).Lat;   % Picked the 3rd region
        idx = find(isnan(x)) ;     % find positions of NaNs
        idx = [1 idx length(x)] ;  % append first and last position
        hold on
        for i = 1:length(idx)-1
            pos = idx(i):idx(i+1) ;
            xi = x(pos) ; yi = y(pos) ;
            % Remove NaN's
            xi(isnan(xi)) = [] ;
            yi(isnan(yi)) = [] ;
            patch(xi,yi,cnty_cdf(j)) ;     % color the area
            c = colorbar;
            c.TickLabelInterpreter = 'latex';
            caxis([0 1.7])
        end
    end
    xlabel('Longitude', 'Interpreter', 'latex'); ylabel('Latitude', 'Interpreter', 'latex')
    title('Failure Rate (FR-2)', 'Interpreter', 'latex', 'FontSize', 16)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)        
    
end