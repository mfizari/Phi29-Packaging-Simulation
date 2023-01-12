function [output_state, transition_time] = StatePath(input_state, rates, subunit)

% This function determines the nucleotide state accessed given the input
% state, the various kinetic rates, and the subunit # (1 is special)
% E->(D,T), D->E, T->(E, Tb), Tb -> Nothing

%input_state = (0,1,2): ADP, Empty, ATP
%Output_state = (0,1,2,3): ADP, Empty, ATP, ATP-tight

%% Define rates and whatnot

k_ADP_release = rates(1);                           %1/s: ADP release
k_ADP_binding = rates(2);                           %1/(s): ADP binding given concentration
k_ATP_binding = rates(3);                           %1/(s): ATP docking given concentration
k_ATP_release = rates(4);                           %1/s: ATP release
k_ATP_tightbind = rates(5);                         %1/s: ATP tight binding
k_ADP_release_sp = rates(6);                        %Special subunit ADP release
k_ATP_binding_sp = rates(7);                      %Special subunit ATP tight binding

%Set rates based on subunit number
if subunit == 1
    k_ADP_release = k_ADP_release_sp;
    k_ATP_binding = k_ATP_binding_sp;  
end


%% Run state selector

if input_state == 0 %if ADP bound, go to empty

    output_state = 1; %Empty
    transition_time = exprnd(1/k_ADP_release);
    
elseif input_state == 1 %if empty, go to ADP bound or ATP bound
    
    t_ED = exprnd(1/k_ADP_binding);
    t_ET = exprnd(1/k_ATP_binding);
    if t_ED > t_ET %ATP binds faster
        output_state = 2;       %ATP bound
        transition_time = t_ET; %ATP binding time
    else
        output_state = 0;       %ADP bound
        transition_time = t_ED; %ADP binding time
    end
    
elseif input_state == 2 %if ATP-bound, lose ATP or tight-bind
    
    t_TE = exprnd(1/k_ATP_release);
    t_TTb = exprnd(1/k_ATP_tightbind);
    if t_TE > t_TTb %ATP tight binds faster
        output_state = 3;        %ATP tightly bound
        transition_time = t_TTb; %ATP binding time
    else
        output_state = 1;       %Empty
        transition_time = t_TE; %ATP unbinding time
    end

end

