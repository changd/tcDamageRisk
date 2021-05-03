clear
dataPath = '/Users/changd/Dropbox (MIT)/2019_DerekAminEmanuelLin/MATLAB_3/matlabData/saturationData/matlabData/';

load([dataPath 'damageRmData.mat']);
totalFailuresRm = totalFailures;

load([dataPath 'damageVmData.mat']);
totalFailuresVm = totalFailures;

load([dataPath 'damageRmVm.mat']);
totalFailuresRmVm = totalFailures;

load([dataPath 'storm.mat']);

Vm = 21:80; Rm = 20:50;
  
vCrit = 20.6; 
S_g = 6.5;

[~,lengthRms] = size(totalFailuresVm); [~,lengthVms] = size(totalFailuresRm);
normalizedDamageVm = (S_g/10)*totalFailuresVm/(prod(size(parameters.latGrid)));
normalizedDamageRm = (S_g/10)*totalFailuresRm/(prod(size(parameters.latGrid)));
normalizedDamageRmVm = (S_g/10)*totalFailuresRmVm/(prod(size(parameters.latGrid)));
normalizedDamageRmVm(1, :) = []; % remove first row to avoid lack of precision for Vm < Vcrit

expected_damages_Vm = zeros(length (Vm), lengthRms);
expected_damages_Rm = zeros(length (Rm), lengthVms);

%% Vm for fixed Rms
figure
subplot(1,2,1)

for i = 1:lengthRms
    for j = 1:length(Vm)
        damage = normalizedDamageVm(j, i);
        
        a = 0;
        b = 0;
    
        for k = 1:(S_g - 1)
            a = a + damage^(k-1) / factorial(k-1);    %left-hand summation of eq.6
        end
        
        for k = 0:(S_g - 1)
            b = b + damage^k / factorial(k);          %right-hand summation of eq.6
        end  
        
        expected_damages_Vm (j, i) = damage * exp(-damage) * a + S_g * (1 - exp(-damage) * b);
    end
    h(i) = plot(Vm, expected_damages_Vm(:, i), 'LineWidth', 2);
    hold on
    
end
    
xlabel('$V_m$','Interpreter','latex'); ylabel('$\bar{\Lambda}_{total}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 20)   
ylim([0 12])
% yticks(0:20:120)

legend([h(1) h(2) h(3) h(4)], {'$R_m$ = 20 km', '$R_m$ = 30 km', ...
     '$R_m$ = 40 km', '$R_m$ = 50 km'}, 'Interpreter', 'latex', 'Location', 'northwest', ...
     'FontSize', 16)


%% Contour plot
[lengthVms, ~] = size(normalizedDamageRmVm); [~,lengthRms] = size(normalizedDamageRmVm); 

Vms = 25:5:80; Rms = 20:2.5:50;


expected_damages = zeros(lengthVms, lengthRms);

for i = 1: lengthRms
    for j = 1: lengthVms
        damage = normalizedDamageRmVm (j, i);
        
        a = 0;
        b = 0;
    
        for k = 1:(S_g - 1)
            a = a + damage^(k-1) / factorial(k-1);    %left-hand summation of eq.6
        end
        
        for k = 0:(S_g - 1)
            b = b + damage^k / factorial(k);          %right-hand summation of eq.6
        end  
        
        expected_damages (j, i) = damage * exp(-damage) * a + S_g * (1 - exp(-damage) * b);
    end
end
        
subplot(1,2,2)
contourf(Vms, Rms, expected_damages')
c = colorbar;
c.TickLabelInterpreter = 'latex';
caxis([0 6.5]);
title('$\bar{\Lambda}_{total}$', 'Interpreter', 'latex')
xlabel('$V_m$', 'Interpreter', 'latex'); ylabel('$R_m$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 20)   
set(gcf,'Renderer', 'painters', 'Position', [10 10 800 300])