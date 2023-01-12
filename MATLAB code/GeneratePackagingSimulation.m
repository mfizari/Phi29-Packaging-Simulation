function [vel, nmin, dwells, bursts, burstsizes] = GeneratePackagingSimulation(rates, ntrials, min_spont)

[dwells, bursts, burstsizes] = Generate_BurstDwells(rates, ntrials, min_spont);


vel = sum(burstsizes)/sum([dwells; bursts]); %vel in bp/s
nmin = mean(dwells)^2 / var(dwells);

mean(dwells)
mean(bursts)
mean(burstsizes)