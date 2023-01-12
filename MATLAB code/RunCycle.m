function [dwell_time, burst_time, burst_size] = RunCycle(rates, k_hydrolysis, k_hydrolysis_spont_sp)

%This function runs 1 mechanochemical cycle and calculates the dwell and
%burst time and the burst size.

%So far, nothing related to filling has been implemented.

% [output_state, transition_time] = StatePath(input_state, rates, subunit)
% The above function determines which path a subunit takes based on its
% current state, and the time. Tightbind is not a valid input state. 
%Input states = 0:2: D,E,T



min_spont_subunits = 3; %Require this many subunits to be ATP-tightbound before spontaneous hydrolysis can occur

Nsubunits = 5;
dwell_time = 0;       %This keeps track of the current cycle dwell time
% cycle_hydrolysis_time = 0;  %This keeps track of the current cycle hydrolysis time
burst_size = 0;             %This keeps track of the current cycle burst size
t_spont_hydrolysis = exprnd(1/k_hydrolysis_spont_sp);
run_cycle = true;

while run_cycle %Running a single cycle of the motor
    
    state = [0,0,0,0,0]; %initialize motor state as all ADP-bound
    for i=1:Nsubunits
        
        run_subcycle = true; %run NT exchange of ith subunit
        
        while run_subcycle
            
            [state(i), transition_time] = StatePath(state(i), rates, i); %generate next state
            
            %Check if spontaneous hydrolysis occurs and ends dwell 
            %If at least the first two have tightly bound ATP and the spont
            %time is shorter than the total dwell time
            if (i > min_spont_subunits) &&  dwell_time+transition_time > t_spont_hydrolysis
                run_cycle = false; %End dwell (outside loop)
                dwell_time = t_spont_hydrolysis;
                burst_size = 2.5*sum(state(2:end) == 3); %each non-special Tb subunit bursts
                run_subcycle = false; %End dwell
            else %Less than 2 subunits have bound ATP or spont hyd. didn't occur fast enough
                
                dwell_time = dwell_time + transition_time;       %keep track of time
                
                if (state(i) == 3) && (i==5) %If cycle is finished
                    burst_size = 10;
                    run_cycle = false;
                    run_subcycle = false;
                elseif (state(i) == 3)      %If given subunit has Tb
                    run_subcycle = false;
                end
                
               
            end
            
            
            
        end
        
    end  
end

%Generate hydryolysis time based on number of subunits firing
burst_time = sum(exprnd(1/k_hydrolysis, 1+burst_size/2.5,1));