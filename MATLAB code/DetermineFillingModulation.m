%% DetermineFillingModulation

%% Experimental data
for z = 1
    
    BurstVsFilling = [76.39247778742734, 13.47384594829316
    78.85957645801322, 21.146369162936764
    81.40924577460083, 32.00013360945957
    83.91161734250784, 38.60444919500301
    86.41278642527891, 44.676999131538565
    89.01997461420268, 64.30022045560828
    91.47665174694367, 50.69744137884963
    94.07181508450797, 65.00300621283989
    96.46175429220389, 71.88723361614005
    99.05751887233615, 86.45868127463429];
    BurstVsFilling(:,2) = BurstVsFilling(:,2)/1e3; %convert to s

    DwellVsFilling = [20.000, 0.106000000
    55.08119073541591, 0.1653471427375004
    65.07979490615375, 0.16346277323357816
    75.09179903780836, 0.1796683509673097
    85.20690842429488, 0.2350660227241005
    92.89327489461488, 0.31166075765612355
    97.88215478815965, 0.5466486139415426];


    VmaxVsFilling = [20, 120.2466168556351
    25, 118.82021004232197
    30, 118.30555410418519
    35, 117.79089816604841
    40, 113.1851193897278
    45, 102.21537175401015
    50, 93.97826427713049
    55, 83.91504258320708
    60, 77.03380531898222
    65, 68.79408537541147
    70, 54.642353309995315
    75, 51.84962641726321
    80, 40.88510371492768
    87.5, 29.915356079210028
    92.5, 19.854746851977694
    97.5, 10.70066356653956];


    KmVsFilling = [20, 26.74666666666666
    25, 28.165517241379305
    30, 27.825287356321834
    35, 28.842298850574707
    40, 27.542068965517238
    45, 24.639999999999997
    50, 22.459770114942526
    55, 21.799540229885054
    60, 21.137471264367814
    65, 18.79632183908046
    70, 15.01609195402299
    75, 16.27402298850575
    80, 15.61379310344828
    87.5, 12.471724137931034
    92.5, 11.010574712643681
    97.5, 9.710344827586209];

    %Renormalize Vmax and Km to account for velocity difference in 2014 vs
    %2009 papers
    VmaxVsFilling(:,2) = VmaxVsFilling(:,2)*(94/120);
    KmVsFilling(:,2) = KmVsFilling(:,2)*(17/26);
    


end

%% Method 1: Determine k_ATP_tightbind from dwell time vs filling

Ntrials = 10e3;

filling = [20;55;65;75;85;92.89; 97.75];
k_ATP_tightbind_array = 1e3*[1,0.070,0.062,0.049,0.0265,0.0151, 0.0036]';
mean_dwell = zeros(length(filling),1);

for i=1:length(filling)

    %Set kinetic rates
    k_hydrolysis = calc_k_hydrolysis(filling(i));
    k_ATP_tightbind = k_ATP_tightbind_array(i);
    
    %Define rate array
    rates = [k_ADP_release, ADP_concentration*k_ADP_binding, ATP_concentration*k_ATP_binding,...
        k_ATP_release,k_ATP_tightbind,k_ADP_release_sp, ATP_concentration*k_ATP_binding_sp ];

    %Preform packaging simulation 
    [~, ~, dwell_durations, ~, ~] = RunPackagingSimulation(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);
%     k_ATP_tightbind_array(i) =  k_ATP_tightbind;
    
%     mean(dwell_durations) - DwellVsFilling(i,2)
    mean_dwell(i) = mean(dwell_durations);
    
end

hold on
plot(filling, mean_dwell, 'ro-')
plot(filling, DwellVsFilling(:,2), 'ko-')
hold off

%% Method 2: Determine k_ATP_tightbind from Vmax/Km data

%Define filling and ATP tightbinding initialization
Ntrials = 5e3;
filling = [20;55;65;75;87.5;92.5; 97.5];
k_ATP_tightbind = zeros(length(filling),1);
k_ATP_tightbind_initial = 1e3*[1,0.070,0.062,0.049,0.0265,0.0151, 0.0036]';

%Reshape experimental data to conform to filling array.
[~, index1, ~] = intersect(VmaxVsFilling(:,1), filling);
[~, index2, ~] = intersect(KmVsFilling(:,1), filling);
Vmax_data = VmaxVsFilling(index1, :);
Km_data = KmVsFilling(index2, :);

ploton = false;

for i=1:length(filling) %for each filling level
    
    tic
    
    %Generate scanning range for ATP tightbind and k_hydrolysis
    midpoint = k_ATP_tightbind_initial(i);
    k_ATP_tightbind_scan = linspace(midpoint-0.25*midpoint, midpoint+0.9*midpoint, 15);
    k_hydrolysis = calc_k_hydrolysis(filling(i));

    
    %Scan and calculate error
    err = zeros(length(k_ATP_tightbind_scan),1);
    for j=1:length(k_ATP_tightbind_scan)
        
        %Define rate array WITHOUT ATP concentration
        rates = [k_ADP_release, ADP_concentration*k_ADP_binding, k_ATP_binding,...
             k_ATP_release,k_ATP_tightbind_scan(j),k_ADP_release_sp, k_ATP_binding_sp];
        
         %Calculate Vmax and Km
        [Vmax_point, Km_point] = CalculateVmaxKm(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);
        
        %Calculate error
        err(j) = ((Vmax_point -Vmax_data(i,2))./(Vmax_data(i,2))).^2 + ((Km_point -Km_data(i,2))./(Km_data(i,2))).^2 ; 
        
        %Plot if desired
        if ploton
           subplot(2,2,1)
           hold on; plot(1,Vmax_point, '.', 'MarkerSize', 15); plot(2,Vmax_data(i,2), '.', 'MarkerSize', 15); hold off
           subplot(2,2,2)
           hold on; plot(1,Km_point, '.', 'MarkerSize', 15); plot(2,Km_data(i,2), '.', 'MarkerSize', 15); hold off
           pause; clf;
        end
        
    end
    
    %Plot 
    plot(k_ATP_tightbind_scan, err, '.', 'MarkerSize', 20)
    pause;
    
    %Assign tightbind value
    [~, index] = min(err);
    k_ATP_tightbind(i) = k_ATP_tightbind_scan(index);
    
    toc

end





%Confirm k_ATP_tightbind by plotting against experimental data
%Generate and plot Vmax, Km, and Vmax/Km to compare to experimental data
for z = 1
    Ntrials = 10e3;
    Vmax_sim = zeros(length(filling),1);
    Km_sim = zeros(length(filling),1);
    mean_dwell = zeros(length(filling),1);

    for i=1:length(filling)

        %Define rate array WITHOUT ATP concentration and k_hydrolysis
        rates = [k_ADP_release, ADP_concentration*k_ADP_binding, k_ATP_binding,...
             k_ATP_release,k_ATP_tightbind(i),k_ADP_release_sp, k_ATP_binding_sp];
        k_hydrolysis = calc_k_hydrolysis(filling(i));

        %Calculate Vmax and Km
        [Vmax_sim(i), Km_sim(i)] = CalculateVmaxKm(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);

        %Calculate mean dwell time WITH ATP_concentration
            ATP_concentration = 500; %uM
            %Define rate array
            rates = [k_ADP_release, ADP_concentration*k_ADP_binding, ATP_concentration*k_ATP_binding,...
                k_ATP_release,k_ATP_tightbind(i),k_ADP_release_sp, ATP_concentration*k_ATP_binding_sp ];

            %Preform packaging simulation 
            [~, ~, dwell_durations, ~, ~] = RunPackagingSimulation(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);
            mean_dwell(i) = mean(dwell_durations);

    end

    %Plot to compare  -pretty close except for the last high-filling Km
    %value??? Also dwells are somehow a bit off....
    subplot(2,2,1)
        hold on
        plot(Vmax_data(:,1), Vmax_data(:,2), '^', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
        plot(filling, Vmax_sim, 'o', 'MarkerFaceColor', 'r' )
        hold off
        xlabel('Filling (%)'); ylabel('V_m_a_x (bp/s)')
    subplot(2,2,2)
        hold on
        plot(Km_data(:,1), Km_data(:,2), 's', 'MarkerFaceColor', 'g', 'MarkerSize', 8)
        plot(filling, Km_sim, 'o', 'MarkerFaceColor', 'r' )
        hold off
        xlabel('Filling (bp/s)'); ylabel('K_m (uM)')
%     subplot(2,2,3)
%         hold on
%         plot(filling, Vmax_sim./Km_sim, 'ro')
%         hold off
    subplot(2,2,4)
        hold on
        bar(filling, DwellVsFilling(:,2))
        plot(filling, mean_dwell, 'ro-')
        hold off
        xlabel('Filling (%)'); ylabel('Dwell time (s)')

end




