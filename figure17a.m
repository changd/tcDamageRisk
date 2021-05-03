clear


%% Determine the best-fit polynomial relating damage to max(V_crit, V_m) and R_m
% damage = c_1 + c_2*R_m*V^c_4 + c_3*R_m^2*V^2c_4
% estimate c_1, c_2, c_3, for fixed c_4 to determine best c_4

load('matlabData/damageRmVm.mat');
totalFailuresRmVm = totalFailures;

load('matlabData/storm.mat');

normalizedDamageRmVm = totalFailuresRmVm/(prod(size(parameters.latGrid)));
normalizedDamageRmVm(1, :) = []; % remove first row to avoid lack of precision for Vm < Vcrit

[lengthVms, ~] = size(normalizedDamageRmVm); [~,lengthRms] = size(normalizedDamageRmVm); 

Vms = 25:5:80; Rms = 20:2.5:50;
c_4s = 1:.01:1.2; c_5s = -.1:.01:0;

degree_error_coeffs = [];
Vs = zeros (size(Vms));
% Velocity matrix with V = max(Vm, vCrit)
 for i = 1:length(Vms)
    V = Vms(i);
        
    Vs(i) = max(V, parameters.vCrit);
 end
 


%% Determine best-fit line

vCrit = parameters.vCrit; lambdaNorm = 3.5*10^-5; nTimes = 121;
degree_error_coeffs = []; rmse = []; cs = [];
pValues = [];

for i = 1:length(c_4s)
    for j = 1:length(c_5s)
        c_4 = c_4s(i); c_5 = c_5s(j);
        data = [];  

        for k = 1:lengthRms
            Rm = Rms(k);

            for l = 1:lengthVms
                V = Vs(l);

                damage = normalizedDamageRmVm (l, k) - lambdaNorm*l*nTimes;

                dataPoints = [Rm*((V-vCrit)/vCrit)^(c_4), Rm^2*((V-vCrit)/vCrit)^(2*c_4), ...
                    -Rm*((V-vCrit)/vCrit)^(c_5), -Rm^2*((V-vCrit)/vCrit)^(2*c_5), damage];
                data = [data; dataPoints];
            end
        end

        mdl = fitlm(data (:, 1:4), data (:, 5));

        coeff = table2array(mdl.Coefficients);
        coeff = coeff (:, 1)'; %[c_1, c_2, c_3]

        degree_error_coeffs = [degree_error_coeffs; coeff];
        pValues = [pValues; mdl.Coefficients.pValue'];
        rmse = [rmse; mdl.RMSE];
        c = [c_4, c_5];
        cs = [cs; c];
    end
end

%%
[~, idx] = min(rmse);
cOptimal = cs(idx,:);
coeffsOptimal = degree_error_coeffs(idx,:);
pValuesOptimal = pValues(idx,:);
vCritVec = vCrit*ones(size(Vs));

figure
idx = 1;
for i = 1:length(Rms)
    Rm = Rms(i);
 
    if mod(i-1,4) == 0
    
        h(idx) = plot([20 Vs], [0; normalizedDamageRmVm(:,i)], 'LineWidth', 4);
        hold on
        plot(Vs, coeffsOptimal(1) + coeffsOptimal(2)*Rm*((Vs-vCritVec)/vCrit).^cOptimal(1) ...
            + coeffsOptimal(3)*Rm^2*((Vs-vCritVec)/vCrit).^(2*cOptimal(1)) ...
            - coeffsOptimal(4)*Rm*((Vs-vCritVec)/vCrit).^cOptimal(2) , '--k', 'LineWidth', 1)
        hold on

        xlabel('$V_m$','Interpreter','latex'); ylabel('$\bar{\Lambda}_{total}$','Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        set(gca, 'FontSize', 15)   
        ylim([0 52])
        idx = idx + 1;
    end
end

set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 20)
legend([h(1) h(2) h(3) h(4)], {'$R_m$ = 20 km', '$R_m$ = 30 km', ...
 '$R_m$ = 40 km', '$R_m$ = 50 km'}, 'Interpreter', 'latex', 'Location', 'northwest')
