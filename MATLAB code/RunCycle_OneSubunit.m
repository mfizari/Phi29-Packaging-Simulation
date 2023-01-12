function [a,b] = RunCycle_OneSubunit(rates, Ntrials)

% Given a rate array, this function runs Ntrials # of simulations to find
% how long it takes for the subunit to reach the tightbind state. 


%% Run cycle to generate tightbind distribution
dist_tightbind = zeros(Ntrials,1);
for i=1:Ntrials

    state = 0;
    run_cycle = true;
    t = 0;

    while run_cycle

        [state, transition_time] = StatePath(state, rates, 2); %generate next state
        t = t + transition_time;

        if state == 3
            dist_tightbind(i) = t;
            run_cycle = false;
        end


    end




end

%%Fit distribution to produce gamma distribution values
pd = fitdist(dist_tightbind,'Gamma');
a = pd.a;
b = pd.b;


