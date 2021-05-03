clear
V_Hol = @(Vm, Rm, B, R) Vm * ((Rm./R).^B .* exp(ones(size(R)) - (Rm./R).^B)).^0.5;

%%
Vms = [25; 37; 46]; Rm = 20:50; R_crits_Rm = zeros(length(Rm), length(Vms));
Vtrans = pos2dist(25,-80,25,-80.1,2); nTimes = 121;
degreesRcrit = zeros(2, length(Vms)); degreesCritZone = zeros(2, length(Vms));
critZoneAreaRm = zeros(length(Rm), length(Vms));

for iter = 1:length(Vms)
    Vm = Vms(iter); B = 1; vCrit = 20.6; R_crit_Rm = zeros(size(Rm));
    R = 0:1500;

    for i = 1:length(Rm)
        V_Hol_Diff = @(R) Vm*((Rm(i)/R)^B * exp(1 - (Rm(i)/R)^B))^0.5 - vCrit;
        R_crit_Rm(i) = fsolve(V_Hol_Diff,51);
        V_HolPlot = V_Hol(Vm, Rm(i), B, R);
    end

    xBasis = [ones(size(Rm')) log(Rm)']; yBasis = log(R_crit_Rm)';
    degree = inv(xBasis'*xBasis)*xBasis'*yBasis;
    
    R_crits_Rm(:,iter) = R_crit_Rm';
    degreesRcrit(:,iter) = degree;
    
    critZoneAreaRm(:,iter) = 2*Vtrans*nTimes*R_crits_Rm(:,iter) + pi*R_crits_Rm(:,iter).^2;    
    xBasis = [ones(size(Rm')) log(Rm)']; yBasis = log(critZoneAreaRm(:,iter));
    degree = inv(xBasis'*xBasis)*xBasis'*yBasis;

    degreesCritZone(:,iter) = degree;
end

%%
Rms = [20; 30; 40; 50]; Vm = 20:80; R_crits_Vm = zeros(length(Vm), length(Rms));
Vtrans = pos2dist(25,-80,25,-80.1,2); nTimes = 121;
degreesRcritVm = zeros(2, length(Rms)); degreesCritZoneVm = zeros(2, length(Rms));
critZoneAreaVm = zeros(length(Vm), length(Rms));

for iter = 1:length(Rms)
    RmIter = Rms(iter); B = 1; vCrit = 20.6; R_crit = zeros(size(Vm));
    R = 0:1500;

    V_Hol_Diff = @(R) Vm(1)*((RmIter/R)^B * exp(1 - (RmIter/R)^B))^0.5 - vCrit;
    R_crit(1) = fsolve(V_Hol_Diff, 60);

    for i = 1:length(Vm)
        V_Hol_Diff = @(R) Vm(i)*((RmIter/R)^B * exp(1 - (RmIter/R)^B))^0.5 - vCrit;
        R_crit(i) = fsolve(V_Hol_Diff,60);
        V_HolPlot = V_Hol(Vm(i), RmIter, B, R);
    end

    xBasis = [ones(size(Vm')) log(Vm)']; yBasis = log(R_crit)';
    degree = inv(xBasis'*xBasis)*xBasis'*yBasis;
    
    R_crits_Vm(:,iter) = R_crit';
    degreesRcritVm(:,iter) = degree;
    
    critZoneArea = 2*Vtrans*nTimes*R_crits_Vm(:,iter) + pi*R_crits_Vm(:,iter).^2;    
    xBasis = [ones(size(Vm')) log(Vm)']; yBasis = log(critZoneArea);
    degree = inv(xBasis'*xBasis)*xBasis'*yBasis;
    
    degreesCritZoneVm(:,iter) = degree;
    critZoneAreaVm(:,iter) = critZoneArea;
end

%%
Vm = 20:80; Rm = 20:50; B = 1; vCrit = 20.6;

Vtrans = pos2dist(25,-80,25,-80.1,2); nTimes = 121;

%% Critical Radius as a function of Vm and Rm
% Initialize matrices
R_crits = zeros(length(Rm), length(Vm)); 

for j = 1:length(Vm)
    for i = 1:length(Rm)
        
        VmPlot = Vm(j); RmPlot = Rm(i);
        
        V_Hol_Diff = @(R) VmPlot *((RmPlot/R)^B * exp(1 - (RmPlot/R)^B))^0.5 - vCrit;  
        R_crits(i, j) = fsolve(V_Hol_Diff,200);  
    end
end

%% Critical Zone Area as a function of Vm and Rm

Area_crits = zeros(length(Rm), length(Vm)); 

for j = 1:length(Vm)
    for i = 1:length(Rm)
        Area_crits(i, j) = 2*Vtrans*nTimes*R_crits(i,j) + pi*R_crits(i,j).^2; 
        
    end
end

%%
figure

for i = 1:length(Rms)
    h(i) = plot(Vm, critZoneAreaVm(:,i), 'LineWidth', 2);
    hold on
    plot(Vm, exp(degreesCritZoneVm(1,i))*Vm.^degreesCritZoneVm(2,i), '--b', 'LineWidth', 2)
    hold on
    xlabel('$V_m$', 'Interpreter', 'latex'); ylabel('$A_{crit}$ (km$^2$)', ...
        'Interpreter', 'latex')
    xlim([min(Vm) max(Vm)])
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 16)
    degreesCritZone(:,i) = degree;
end
legend([h(1), h(2), h(3), h(4)], {'$R_m$ = 20 km', '$R_m$ = 30 km', ...
    '$R_m$ = 40 km', '$R_m$ = 50 km'}, 'Location', 'northwest', 'Interpreter', 'latex')

set(gcf,'Renderer', 'painters', 'Position', [10 10 400 300])