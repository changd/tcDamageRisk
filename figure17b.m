clear

%% Determine the best-fit polynomial relating damage to max(V_crit, V_m) and R_m

load('matlabData/LossRmVm.mat');
totalLossRmVm = totalCost;

load('matlabData/storm.mat');

normalizedLossRmVm = totalLossRmVm/(prod(size(parameters.latGrid)));
[lengthVms, lengthRms] = size(normalizedLossRmVm); 

Vms = 20:5:80; Rms = 20:5:50; Vms(1) = 20.6;
d6s = 1.2:.01:2; 
vCrit = 20.6;

%% Determine best-fit line
degree_error_coeffs = []; rmse = []; ds = []; pValues = [];

for i = 1:length(d6s)
    
    d6 = d6s(i); 
    data = [];  
    losses = [];
    for k = 1:lengthRms
        Rm = Rms(k);

        for l = 1:lengthVms
            Vm = Vms(l);
            gVm = (Vm-vCrit)/vCrit;

            loss = normalizedLossRmVm (l, k);

            dataPoints = [Rm * gVm^(d6), Rm^2 * gVm^(2*d6), Rm^3 * gVm^(3*d6), Rm^4 * gVm^(4*d6), ...
                Rm^2 * gVm^(d6), Rm^3 * gVm^(d6), Rm^3 * gVm^(2*d6), Rm^4 * gVm^(2*d6), ...
                Rm, Rm^2, Rm^3, Rm^4];
            dataPoints = real(dataPoints);

            data = [data; dataPoints];
            losses = [losses; loss];
        end
    end

    mdl = fitlm(data, losses);
    coeff = table2array(mdl.Coefficients);
    coeff = coeff (:, 1)'; 

    degree_error_coeffs = [degree_error_coeffs; coeff];
    pValues = [pValues; mdl.Coefficients.pValue'];
    rmse = [rmse; mdl.RMSE];
    ds = [ds; d6];
end

%%
[~, idx] = min(rmse);
dOptimal = ds(idx,:);
coeffsOptimal = degree_error_coeffs(idx,:);
pValuesOptimal = pValues(idx,:);
vCritVec = vCrit*ones(size(Vms));

%%
figure
idx = 1;
gVms = (Vms-vCritVec)/vCrit;
onesv = ones(size(vCritVec));

for i = 1:lengthRms
    Rm = Rms(i);

    if mod(i-1,2) == 0
        Vms(1) = 20;
        h(idx) = plot(Vms, normalizedLossRmVm(:,i), 'LineWidth', 3); 
        hold on
        
        plot(Vms, coeffsOptimal(1) + ...
            coeffsOptimal(2)*Rm*gVms.^dOptimal(1) + coeffsOptimal(3)*Rm^2*gVms.^(2*dOptimal(1)) + ...
            coeffsOptimal(4)*Rm^3*gVms.^(3*dOptimal(1)) + coeffsOptimal(5)*Rm^4*gVms.^(4*dOptimal(1)) + ...
            coeffsOptimal(6)*Rm^2*gVms.^dOptimal(1) + coeffsOptimal(7)*Rm^3*gVms.^(dOptimal(1)) + ... 
            coeffsOptimal(8)*Rm^3*gVms.^(2*dOptimal(1)) + coeffsOptimal(9)*Rm^4*gVms.^(2*dOptimal(1)) + ...
            coeffsOptimal(10)*Rm*onesv + coeffsOptimal(11)*Rm^2*onesv + ...
            coeffsOptimal(12)*Rm^3*onesv + coeffsOptimal(13)*Rm^4*onesv, '--w', 'LineWidth', 1);
        hold on

        xlabel('$V_m$','Interpreter','latex'); ylabel('$\bar{L}_{total}$','Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        set(gca, 'FontSize', 15)   
        ylim ([0 2500])
        idx = idx + 1;
    end

end

set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 20)
legend([h(1) h(2) h(3) h(4)], {'$R_m$ = 20 km', '$R_m$ = 30 km', ...
 '$R_m$ = 40 km', '$R_m$ = 50 km'}, 'Interpreter', 'latex', 'Location', 'northwest')
