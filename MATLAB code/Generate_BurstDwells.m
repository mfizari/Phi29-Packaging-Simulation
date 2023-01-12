function [dwells, bursts, burstsizes] = Generate_BurstDwells(rates, ntrials,min_subunits_for_spont)


dwells = zeros(ntrials,1);
burstsizes = zeros(ntrials,1);
bursts = zeros(ntrials,1);

for k=1:ntrials

    %Dwell for special subunit
    dwell_time = Generate_SingleSubunit_Tightbindtime(rates,1) ;

    %Dwell for next two subunits
    dwell_time = dwell_time + Generate_SingleSubunit_Tightbindtime(rates, min_subunits_for_spont-1);

    %Calculate spont hyd. time
    t_spont = exprnd(1/rates.k_sponthydrolysis);

    %Dwell for remaining subunits
    n_subunits_firing=0;
    for i=min_subunits_for_spont+1:5
        t_dwell_subunit=Generate_SingleSubunit_Tightbindtime(rates,1);
        if t_spont < t_dwell_subunit
            n_subunits_firing=i-2;
            continue;
        else
            dwell_time = dwell_time + t_dwell_subunit;
        end
    end

    %Record hydrolysis time, dwell time, and number of subunits firing/
    if n_subunits_firing==0 %spont hyd didnt occur
        n_subunits_firing = 4; %all 4 fired
        dwells(k) = dwell_time;
        burstsizes(k) = 2.5*n_subunits_firing;
        bursts(k) = exprnd(1/rates.k_hydr);
    else
        dwells(k) = dwell_time;
        burstsizes(k) = 2.5*n_subunits_firing;
        bursts(k) = (n_subunits_firing/5)*exprnd(1/rates.k_hydr);

    end
end


%% Define functions
function transition_time = GenerateEmptyToTightbindTime(rates)

p = (1/rates.k_ATP_unbind)/((1/rates.k_ATP_unbind)+(1/rates.k_ATP_tightbind)); %probability to go from T->Tb
n_steps=1+geornd(p); %random sample from dist - number of trials req. to get to Tb from T. Each trial 

%Number of E->T is n_steps
sum_tET = sum(exprnd(1/rates.k_ATP_bind, n_steps,1));
%Number of T->E is n_steps-1
sum_tTE = sum(exprnd(1/rates.k_ATP_unbind, n_steps-1,1));
%Number of T->Tb is 1
sum_tTTb = exprnd(1/rates.k_ATP_tightbind);

transition_time = sum_tET+sum_tTE+sum_tTTb;

end

function transition_time = Generate_SingleSubunit_Tightbindtime(rates, nsubunits) 

transition_time=0;
for j=1:nsubunits
    transition_time =  transition_time + exprnd(1/rates.k_ADP_unbind) +...
        GenerateEmptyToTightbindTime(rates);
end

end


end
