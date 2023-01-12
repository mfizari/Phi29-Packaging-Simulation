%%PackagingSimulation_AddNoise

cfold = 'E:\Phi29 motor simulation\Code';

%% Define kinetic rates

%Define filling and filling dependent kinetic rates
filling = 90;
k_hydrolysis = calc_k_hydrolysis(filling);          %1/s: hydrolysis rate
k_ATP_tightbind = calc_k_ATP_tightbind(filling);    %1/s: ATP tight binding

%Define nucleotide concentrations
ATP_concentration = 10;                    %uM
ADP_concentraion  = 1;                      %uM

%Assemble rate array
rates = [k_ADP_release, ADP_concentraion*k_ADP_binding, ATP_concentration*k_ATP_binding,...
        k_ATP_release,k_ATP_tightbind, ADP_concentraion*k_ADP_release_sp, ATP_concentration*k_ATP_binding_sp ];
    
%This occurs when at least 2 subunits (special + next one) are ATP bound.
k_hydrolysis_spont_sp = 0.4*4;      %1/s: Special subunit spontaneous hydrolysis 


%% Run packaging simulation

%Data aquisition parameters
sr = 1e3;           %Sample rate for recording 
isr = 1/sr;         %Time step for recording data
Ntrials = 0.5e3;      %Number of dwell/burst cycles. Len.pack = 10*Ntrials (bp)

%Initialize values
dwell_durations = zeros(Ntrials,1);
burst_durations = zeros(Ntrials,1);
burst_sizes = zeros(Ntrials,1); 
tl = [0];      %#ok<NBRAK> %tether length array - bp

for i=1:Ntrials

    % Generate burst and dwell duration for time step
    [dwell_durations(i), burst_durations(i), burst_sizes(i)] = RunCycle(rates, k_hydrolysis, k_hydrolysis_spont_sp);

    % Round dwell and burst to match timestep value
    dwell_durations(i) = isr*floor(dwell_durations(i)/isr);
    burst_durations(i) = isr*floor(burst_durations(i)/isr);

    % Number of time steps in each dwell and burst
    nstep_dwell = floor(dwell_durations(i)/isr);
    nstep_burst = floor(burst_durations(i)/isr);

    % Generate tether lengths in dwell and burst and append
    tl_dwell = tl(end)*ones(nstep_dwell,1);
    tl_burst = linspace(tl(end), tl(end)+burst_sizes(i), nstep_burst)';
    tl = [tl; tl_dwell; tl_burst];


end

%% Add experimental noise to the packaging trace

fname = 'Noise_pd_Win0.1_Slide1.mat'; %this file contains a distribution fitted to the noise of a 10kb control tether with a given window size/slidefactor
win = str2double(fname(13:15));
units = 'bp';
tl_noise = AddNoiseToSimTrace(cfold, tl, fname, win,units);

t = (1:length(tl_noise))/1e3; 

hold on
plot(t,smoothdata(tl_noise, 'gaussian', 1)); 
plot(t, tl, 'r.')
xlabel('Time (s)'); ylabel('Length packaged (bp)')



%% Compare to a typical packaging trace

compare_on = false;

if compare_on
    %Load packaging trace from ejection experiment
    dfold = 'F:\phi29 ejection\Selected high filling\High fill - high force\5 pN alldata\ej - selected';
    load('13-12-57.mat', 'lpkg');
    
    %Convert lpkg to length_packaged (positive values) and only first 1kbp
    lpkg = -1e3*(lpkg-mean(lpkg(1:100)));
    lpkg = lpkg(lpkg <= 1e3);
    
    
    %Polyfit both and subtract for noise only
    t1 = (1:length(lpkg))'; t2 = (1:length(tl_noise))';
    f1 = polyfit(t1, lpkg,1);
    f2 = polyfit(t2, tl_noise,1);
    f1 = t1*f1(1) + f1(2);
    f2 = t2*f2(1) + f2(2);
    
    hold on
    plot(t1/1e3, lpkg-f1)
    plot(t2/1e3, tl_noise-f2)
    hold off
    title('Comparison of expt. packaging noise to simulation+noise')
    xlabel('Time (s)'); ylabel('Tether length - poly. fit');
    set(gca,'FontSize',12)
end

%% 







