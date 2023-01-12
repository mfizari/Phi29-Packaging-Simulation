function err = ProduceATPFitErr(input_rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials, ATP_concentration_array,ADP_concentration, experimental_data, ploton)


%This function calculates the squared error when comparing the mean
%dwell([ATP]) and the nmin([ATP]) from experiments to that calculated from
%the motor simulation.

%NOTE: [ATP] in experimental data must be rounded
%(5,10,25,50,100,250,500,1000)

%input_rates = [k_ADP_release, k_ADP_binding, k_ATP_binding,...
%               k_ATP_release,k_ATP_tightbind, k_ADP_release_sp, k_ATP_binding_sp];
%Note: input_rates does NOT include ATP concentration or ADP concentration.

% k_hydrolysis = motor hydroylsis rate
% k_hydrolysis_sp = spontaneous hydrolysis rate
% Ntrials = # of packaging cycles to run. Length packaged = ~10*Ntrials
% ATP_concentration_array = [ATP] to scan
% ADP_concentration = given ADP concentration
% experimental_data = cell array:
%    experimental_data{1} = mean dwell time vs [ATP] (nx2 matrix)
%    experimental_data{2} = nmin vs [ATP]            (nx2 matrix)

%% Re-size experimental data to match the ATP_concentration_array

[~, index1, ~] = intersect(experimental_data{1}(:,1), ATP_concentration_array);
[~, index2, ~] = intersect(experimental_data{2}(:,1), ATP_concentration_array);

experimental_data{1} = experimental_data{1}(index1, :);
experimental_data{2} = experimental_data{2}(index2, :);


%% Initialize observables
Nmin = zeros(length(ATP_concentration_array),1);
mean_dwell = zeros(length(ATP_concentration_array),1);

%Multiply ADP binding rate constants by ADP_concentration
input_rates(2) = ADP_concentration*input_rates(2);
input_rates(5) = ADP_concentration*input_rates(5);

%Calculate values for each ATP concentration
for k=1:length(ATP_concentration_array)
    
    %Multiply ADP binding rate constants by ADP_concentration
    input_rates(3) = ATP_concentration_array(k)*input_rates(3);
    input_rates(7) = ATP_concentration_array(k)*input_rates(7);
    
    
    [~, nmin, dwell_durations, ~, ~] = RunPackagingSimulation(input_rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);
    
    
    mean_dwell(k) = mean(dwell_durations);
    Nmin(k) =  nmin;
    
end


%% Calculate the error in the two values

err1 = (((experimental_data{1}(:,2) - mean_dwell))./(experimental_data{1}(:,2))).^2;
err2 = (((experimental_data{2}(:,2) - Nmin))./(experimental_data{2}(:,2))).^2;

err = 0.9*sum(err1) + 0.01*sum(err2);

if ploton
    

    subplot(2,2,1) %compare mean dwell vs [ATP]
        hold on
        plot(experimental_data{1}(:,1), experimental_data{1}(:,2)/1e3, 'ko-')
        plot(ATP_concentration_array, mean_dwell, 'ro-')
        hold off
        set(gca,'XScale', 'log')
    subplot(2,2,2) %compare nmin vs [ATP]
        hold on
        plot(experimental_data{2}(:,1), experimental_data{2}(:,2), 'ko-')
        plot(ATP_concentration_array, Nmin, 'ro-')
        hold off
        set(gca,'XScale', 'log')
    
        pause;
end


