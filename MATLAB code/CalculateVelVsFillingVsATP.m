%% Analysis_CalculateVelVsFillingVsATP

%Using the motor simulation I have and the parameters that most closely
%match up to Bustamante group data, this script calculates velocity vs
%filling curves for different ATP concentrations and compares to published
%data from Liu et al 2014



%% Run packaging simulation

filling = [30, 40:5:90]';
veldata = zeros(length(filling), length(ATP_concentration));

for i=1:length(ATP_concentration)
    
    for j=1:length(filling)
        
        %Set kinetic rates
        k_hydrolysis = calc_k_hydrolysis(filling(j));
        k_ATP_tightbind = calc_k_ATP_tightbind(filling(j));
        
         %Define rate array
        rates = [k_ADP_release, ADP_concentration*k_ADP_binding, ATP_concentration(i)*k_ATP_binding,...
        k_ATP_release,k_ATP_tightbind,k_ADP_release_sp, ATP_concentration(i)*k_ATP_binding_sp ];
        
        %Perform packaging simulation 
        [v, ~, ~ , ~, ~] = RunPackagingSimulation(rates, k_hydrolysis, k_hydrolysis_spont_sp, Ntrials);
        
        %Add to data matrix
        veldata(j,i) = v;
    end
        

end


%% Load curves from Liu et all 2014 into veldata_paper matrix

cfold = '';
dfold = '';
cd(dfold); files = dir('*csv'); filenames = {files.name};
filling_paper = [30, 40:5:90]';
ATP_concentration_paper = flip(ATP_concentration);
veldata_paper = zeros(length(filling_paper), length(ATP_concentration_paper));
for k=1:length(filenames)
    
    x = load(filenames{k});
    x = x(:,2); x = flip(x);
    
    
    veldata_paper(k, :) = x';
    
    
    
end

%Motor velocity vs filling from 2015 Berndsen et al
data = [3.913008968965798, 106.49416620366696
6.884673386776736, 105.83062147789508
11.338995158345906, 104.53686800539725
16.45368680053973, 105.3178823716168
20.885784586078262, 101.93507421223907
26.156044130486553, 97.33947138661799
31.273910627827615, 98.41892213667751
36.76323517739503, 94.41542979601556
41.37947456147313, 88.34193189935709
46.04968648305422, 87.3418525279784
51.31994602746249, 82.74624970235732
56.14890070640527, 76.66798952297802
60.97785538534805, 70.58972934359869
65.53059766648147, 58.54750377014048
70.33097864909912, 49.78331613620128
75.59171362806572, 44.292404159060254
80.37622033494722, 34.0360346059211
85.40519088816572, 26.759266608460976
90.42463687594253, 18.58718945948091
94.61227081514404, 12.22477974442414
98.84435272640684, 10.040479403127222
3.9257083895547353, 107.6879117390269
6.681482657353762, 106.73069291213588
11.358044289229309, 106.32748630843717
16.67910151599334, 106.5068656242559
21.111199301531872, 103.12405746487816
26.381458845940156, 98.52845463925709
31.492975632986745, 99.01103262163663
36.56956901341377, 96.21081038177633
41.39217398206206, 89.53567743471704
46.06238590364315, 88.53559806333836
50.69767441860466, 84.25271846971981
56.16477498214144, 78.16017144217795
60.777839511072315, 71.78823716167949
65.75601238193507, 59.736487022779585
70.77228351456465, 51.26597348995952
75.39487260893722, 45.78934836098102
80.39209461068339, 35.528216525121024
85.41471545360743, 27.654575759980943
90.22144614651955, 19.48726089372171
94.40908008572109, 13.124851178664969
98.86340185729024, 11.831097706167142];

fill = data(1:21, 1);
v = data(1:21,2); verr = data(22:end, 2)-v;
alpha = mean(v(fill>30 & fill < 35));
v = v/alpha; verr = verr/alpha;




data_carlos = [9.377755856833822, 0.5313253012048192
24.331914632168576, 0.26385542168674697
48.41405562237217, 0.18674698795180722
98.32054138576147, 0.14096385542168677
247.69229211342363, 0.11686746987951802
500.9677459532117, 0.10481927710843364
1012.6360692146282, 0.09759036144578304];

data_expt = [9.5938102016134, 0.3433734939759036
23.84654643531229, 0.23012048192771084
48.58400824178851, 0.1578313253012048
96.07828400393093, 0.13132530120481922
241.90220241737018, 0.11204819277108424
489.5428838583685, 0.0951807228915662
964.7188422032918, 0.09759036144578304];


data_sim = [9.849296233331778, 0.32650602409638557
23.2755083854183, 0.23012048192771084
48.54146423720192, 0.1650602409638554
100.73272372299547, 0.14096385542168677
247.69229211342363, 0.11686746987951802
488.82862157441394, 0.10722891566265058
989.253387222846, 0.09036144578313243];

ATP = [10,25,50,100,250,500,1e3]';

t1 = data_expt(:,2);
t2 = data_sim(:,2);
t3 = data_carlos(:,2);







%% Plot velocity vs filling curves

FS = 13;
LW = 1.1;
MS = 4;


subplot(2,2,1) %Plot the other stuff (dwell time simulation)

    hold on
    plot(ATP, t1, 'ko-', 'LineWidth', LW)
    plot(ATP, t2, 'o--', 'LineWidth', LW)
    plot(ATP, t3, 'o--', 'LineWidth', LW, 'Color', [0.4940, 0.1840, 0.5560]	)
    hold off
    set(gca,'XScale', 'log')
    xlabel('[ATP] ({\mu}M)')
    ylabel('Dwell time (s)')
    set(gca,'FontSize', FS); box on;
    legend({'Experiment: Moffit et al. 2009','Current simulation','Simulation: Chistol et al. 2012'})

subplot(2,2,2) %Plot the velocity vs filling curves for different ATP concentrations!

    % Plot expt data
    select = fill >= 30 & fill <= 90;
    hold on
        errorbar(fill(select), v(select), verr(select), 'ks-', 'MarkerSize', MS,...
            'LineWidth', LW)
    hold off

    % Plot Velocity vs filling curves from simulation
    for i=1:size(veldata,2)
        hold on
        plot(filling, veldata(:,i)./max(veldata(:,i)),'.--')
        hold off
        pause;
    end
    
    
    
    box on;
    xlabel('Capsid filling (%)');
    ylabel('Normalized motor velocity')
    set(gca,'FontSize', FS)




    



