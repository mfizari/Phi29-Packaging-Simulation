%% Phi29 Monte Carlo

%Rates are taken from Row S4J in the SI of 
%"High Degree of Subunit Coordination (2012)

%Spontaneous hydrolysis rates are taken from the 2018 PNAS paper

%% Define kinetic rates

rates = [k_ADP_release, ADP_concentraion*k_ADP_binding, ATP_concentration*k_ATP_binding,...
        k_ATP_release,k_ATP_tightbind, ADP_concentraion*k_ADP_release_sp, ATP_concentration*k_ATP_binding_sp ];

%Hydrolysis rate PER subunit: Should depend on filling!
k_hydrolysis = 5/0.01;
    
%This occurs when at least 2 subunits (special + next one) are ATP bound.
%Burst size ranges from 2.5bp to 10bp
k_hydrolysis_spont_sp = 4;      %1/s: Special subunit spontaneous hydrolysis 

% What else to include?
%1. Tight binding modulation with filling (try to reproduce data)
%2. Modulation of spontaneous hydrolysis rate that changes the average
%   burst size as a function of filling??? Or does it only change due to the 
%   change in rotation?(see Liu 2014)
%3. Modulation of the regular hydrolysis rate due to increasing F_int
%4. Presence of stochastic LLPs (characterized in Liu 2014)

%Question: Is the spontaneous hydrolysis rate measured from the ADP-bound 
%state or only once the first two subunits bind ATP??



%% Running packaging simulation

Ntrials = 1e4;                              %Number of burst/dwell cycles to run
dwell_durations = zeros(Ntrials,1);
burst_durations = zeros(Ntrials,1);
burst_sizes = zeros(Ntrials,1); 

for i=1:Ntrials
    [dwell_durations(i), burst_durations(i), burst_sizes(i)] = RunCycle(rates, k_hydrolysis, k_hydrolysis_spont_sp);
end

tpkg = cumsum(dwell_durations);
lpkg = cumsum(burst_sizes);


%Plot stuff
subplot(2,2,1)
    plot(tpkg, lpkg);
    title(num2str(range(lpkg)/range(tpkg)))
    xlabel('Time (s)'); ylabel('DNA packaged (bp)')
subplot(2,2,2)
    histogram(burst_sizes)
    xlabel('Burst size (bp)'); ylabel('Count')
subplot(2,2,3)
    histogram(1e3*dwell_durations)
    xlabel('Dwell duration (ms)'); ylabel('Count')
subplot(2,2,4)
    histogram(1e3*burst_durations)
    xlabel('Burst duration (ms)'); ylabel('Count')