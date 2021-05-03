clear %vars

cdf_times = [33 34 38 61 64 88];  %these are the labels for the times in the cdf files
outage_times = ["1010_1540", "1010_1635", "1010_1950", "1011_1940", "1011_2200", "1012_2305"];  %times in the power outage files
real_times = ["October 10, 16:35", "October 10, 19:50", "October 11, 19:40", "October 11, 22:00", "10/12 11:05 PM"];

alphabet = 'ABCDEF';  %there should to be a better way to say what excel column I want 

dataPathStem = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/MATLAB_3/outageData/';

pValues = zeros(length(outage_times), 1);
% for i = 1:2  %do once with cdf, and once with cumulative velocity exceedance on the x-axis
for i = 2  %do once with cdf, and once with cumulative velocity exceedance on the x-axis
    if i == 1
        dataPath = [dataPathStem, 'regression_data/cdf/'];
    else
        dataPath = [dataPathStem, 'regression_data/cumulative_velocity/'];
    end
    
    for time_step = 2:length(outage_times)

        %raw power outage for each county
        n = readmatrix([dataPathStem, 'Guam Data CSV Files/excel/Power_Outages_2018', convertStringsToChars(outage_times(time_step))], 'Range', 'A2:L38');

        %number of housholds without power; this is the observed response
        number_out = n(:,5);   %considered "number of successes" in nb
        pct_out = n(:,6);      %used in the scatterplot

        %these are the predictors
        if i == 1
            % Get average county failure rate
            cnty_cdf_avg = readmatrix([dataPathStem, 'County_CDF_Calculations/county_cdf_t=' num2str(cdf_times(time_step)) '.xlsx'])';  %average cdf across the counties' grids
            x1 = cnty_cdf_avg;
        else
            % Get cumulative velocity metric
            cumul_vel = readmatrix([dataPathStem, 'County_Velocity_Calculations/county_velocity_cumulative.xlsx'])';
%             cumul_vel = readmatrix(['\Users\HP\Dropbox (MIT)\urop2020_sharedMaterialsDerek\County_Velocity_Calculations\county_velocity_cumulative.xlsx'])';
            x1 = cumul_vel(:, cdf_times(time_step));
        end
        
        % Get other predictors (?) (county area + number of households)
        cnty_area = readmatrix([dataPathStem, 'FLCountyLatLonArea_Excel.xlsx'], 'Range', 'D2:D38');
        total_households = n(:,4);   %considered "number of trials" in nb

        %sort the rows by increasing cdf so the line plot doesn't look funky
        [x1, idx] = sortrows(x1);
        number_out = number_out(idx,:);
        pct_out = pct_out(idx,:);
        cnty_area = cnty_area(idx,:);
        total_households = total_households(idx,:);


         %% binomial distribution, using logit link
         
         if time_step == 2 && i == 1
            idx = find(x1 < 0.4);
            x1 = x1(idx);
            cnty_area = cnty_area(idx);
            number_out = number_out(idx);
            total_households = total_households(idx);
            pct_out = pct_out(idx);
        end
        
         % input
         X = [x1, cnty_area];
         X = x1;
         
         % glmfit using binomial distribution; dispersion parameter
        [b, dev, stats] = glmfit(X, [number_out, total_households], 'binomial', 'estdisp', 'on');
        pValues(time_step) = stats.p(2);
        
        % gets model-predicted outages
        yfit = glmval(b(1:2,:),X(:,1), 'logit');  %this plots yfit automatically
        
        % write model output
        writematrix([outage_times(time_step); yfit], [dataPath 'binomial_reg_summary.xlsx'], 'Sheet', ...
            'Estimated Outage %', 'Range', [alphabet(time_step) ':' alphabet(time_step)]);


        %% Way 2, with binomial distribution
        
        [phat,pci] = binofit(round(yfit .* total_households), total_households);
        writematrix([outage_times(time_step); pci(:,2) - pci(:,1)], [dataPath 'binomial_reg_summary.xlsx'], ...
            'Sheet', 'Error Bar Magnitudes', 'Range', [alphabet(time_step) ':' alphabet(time_step)]);
        

        %% Make Scatterplot + Logistic Regression Curve
        figure
        plot(x1, 100*pct_out, 's', 'MarkerSize', 9, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6])
        hold on
        errorbar(x1, 100*phat, 100*(phat - pci(:,1)), 100*(pci(:,2) - phat), '-r', 'LineWidth', 1.5);
        
        if i == 1
            xlabel('Failure Rate (FR-2)', 'Interpreter', 'latex'); 
        else
            xlabel('Cumulative Velocity', 'Interpreter', 'latex');
        end
        ylabel('Outages per 100 Households', 'Interpreter', 'latex');
        
        legend({'Observed', 'Estimated'}, 'Interpreter', 'latex', 'location', 'northwest');
        
        set(gca,'TickLabelInterpreter','latex')
        set(gca, 'FontSize', 16)

        set(gcf,'Renderer', 'painters', 'Position', [10 10 500 375])
        
        
    end
    
end