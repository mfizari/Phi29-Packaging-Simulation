function [v, nmin, dwell_durations, burst_durations, burst_sizes] = RunPackagingSimulation(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials)

dwell_durations = zeros(Ntrials,1);
burst_durations = zeros(Ntrials,1);
burst_sizes = zeros(Ntrials,1); 

for i=1:Ntrials
    [dwell_durations(i), burst_durations(i), burst_sizes(i)] = RunCycle(rates, k_hydrolysis, k_hydrolysis_spont_sp);
end

cycle_times = dwell_durations+burst_durations;
v = sum(burst_sizes)/range(cumsum(cycle_times)); %bp/s
nmin = mean(dwell_durations)^2 / var(dwell_durations);