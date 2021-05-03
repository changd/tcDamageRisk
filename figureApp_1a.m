clear

%% Determine the best-fit polynomial relating R_crit to max(V_crit, V_m) for given R_m values
% R_crit = = A * R_m^C * [max(V_crit, V_m)]^B     find appropriate A&B&C

% Majority of code is from exampleCritZoneSizeVm.m


V_Hol = @(Vm, Rm, B, R) Vm * ((Rm./R).^B .* exp(ones(size(R)) - (Rm./R).^B)).^0.5; 

%% Plot critical radius as a function of Vm
% Plot four curves -- corresponds to four different values of Rm

Rms = [20; 30; 40; 50]; Vm = 21:80; 

% Initialize matrices
R_crits = zeros(length(Vm), length(Rms)); degreesRcrit = zeros(2, length(Rms)); Vs = zeros(size(Vm));
R_crits_inner = zeros(length(Vm), length(Rms)); 

% Additional parameters
Vtrans = pos2dist(25,-80,25,-80.1,2); nTimes = 121; vCrit = 20.6; 

% Fill velocity matrix with max (Vm, vCrit)
for k = 1:length(Vm)
    V = Vm(k);
    Vs(k) = max(V, vCrit);
end
 
% Initialize data matrix [Vs, Rs, R_crit]
data = []; 

%% Calculate R_crit values for all parameters
for iter = 1:length(Rms)
    
    % Set Holland parameters
    Rm = Rms(iter); B = 1; 

    % Calculate critical radius at all values of Vm
    R = 0:1500;
    R_crit = zeros(size(Vm));
    R_crit_inner = zeros(size(Vm));
    
    for i = 1:length(Vm)
        V_Hol_Diff = @(R) Vm(i)*((Rm/R)^B * exp(1 - (Rm/R)^B))^0.5 - vCrit;  
        
        R_crit(i) = fsolve(V_Hol_Diff,51);
        R_crit_inner(i) = fsolve(V_Hol_Diff,Rm-5);
        dataPoints = [Vm(i), Rms(iter), R_crit(i), R_crit_inner(i)];
        data = [data; dataPoints];
    end

    R_crits(:,iter) = R_crit'; R_crits_inner(:,iter) = R_crit_inner';
end

%% Calculate best-fit line
X = data (:, 1:2);
Y = data (:, 3);
mdl = fitlm(log(X), log(Y));

degree = table2array(mdl.Coefficients);         %[coefficient, degree of V, degree of R]
degree = degree (:, 1)';
degree(1) = exp(degree(1));
degrees = degree;

%% Plot Lines

for iter = 1:length (Rms)

    % Plot critical radius as a function of Vm
    h(iter) = plot(Vm, R_crits (:, iter), 'LineWidth', 2);
    hold on
    
    % Plot polynomial fit
%     plot(Vm, degree(1)*(Vs-vCrit).^2*Rms(iter), '--b', 'LineWidth', 2);
    plot(Vm, degree(1)*Vs.^degree(2)*Rms(iter)^degree(3), '--b', 'LineWidth', 2);
    hold on
    
    %title ('R_{crit} = A*R_m*max(V_{crit}, V_m)^B')
    xlabel('$V_m$', 'Interpreter', 'latex'); ylabel('$R_{crit}$', 'Interpreter', 'latex')
    xlim([min(Vm) max(Vm)])
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    
end

legend([h(1), h(2), h(3), h(4)], {'$R_m$ = 20 km', '$R_m$ = 30 km', ...
    '$R_m$ = 40 km', '$R_m$ = 50 km'}, 'Location', 'northwest', 'Interpreter', 'latex')
