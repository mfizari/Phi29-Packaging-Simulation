function [Vmax, Km] = CalculateVmaxKm(rates, ntrials, min_spont)

%Initialize 
ATP_concentration = [10,25,50,100,250,500,1000,1500]; %uM
vel = zeros(length(ATP_concentration),1);
initial_k_ATP_bind = rates.k_ATP_bind;


%Run packaging simulation for each ATP concentration!
for k=1:length(ATP_concentration)
   
    %New ATP binding rate
    rates.k_ATP_bind = ATP_concentration(k)*initial_k_ATP_bind;

    [v, ~, dwells, ~, ~] = GeneratePackagingSimulation(rates, ntrials, min_spont);
    vel(k)=v;
end

%% Fit curve to get Vmax and Km values

MM_fitfun = @(x,xdata)x(1)*xdata./(x(2) + xdata);
x0 = [100,20];
MM_fit_parameters = lsqcurvefit(MM_fitfun,x0,ATP_concentration',vel);

Vmax = MM_fit_parameters(1);
Km = MM_fit_parameters(2);

