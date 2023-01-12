function [fitvals, scan_array] = ScanRate(rates,ATP_concentration, ADP_concentration, rate_index, scan_bounds, Nscan, Ntrials)

%rates is given WITHOUT nucleotide concentrations!

%You give the standard rate array, and ATP/ADP concentration
%You specify, using rate_index, which value you are scanning
%Code generates a scan array for the value you are scanning, log spaced
%For each scan value, a new rate array is constructed and dist is fit!

%rate_index: [1:7]: 
% 1: k_ADP_release
% 2: k_ATP_binding
% 3: k_ATP_release
% 4: k_ATP_tightbind
% 5: k_ADP_release_sp
% 6: k_ATP_binding_sp
% 7: ATP_concentration

%Shift rate_index so that it corresponds to the correct thingy
if rate_index > 1
    rate_index = rate_index + 1;
end


%Generate scan_array
lb = scan_bounds(1); ub = scan_bounds(2);
scan_array = logspace(log10(lb),log10(ub),Nscan)';

%Preform scan!
rates_scan = rates; 
rates_scan(2) = ADP_concentration*rates(2); %dont ever scan ADP binding stuff
fitvals = zeros(length(scan_array),2);

for i=1:length(scan_array)

    if rate_index == 8 %You are scanning ATP_concentration, k_ATP_b fixed

        rates_scan(3) = rates(3)*scan_array(i); %regular 
        rates_scan(7) = rates(7)*scan_array(i); %special subunit

    elseif rate_index == 3 %You are scanning ATP binding rate, [ATP] fixed

        rates_scan(3) = scan_array(i)*ATP_concentration; 
        rates_scan(7) = rates(7)*ATP_concentration;

    elseif rate_index == 7 %You are scanning sp ATP binding, [ATP] fixed

        rates_scan(3) = rates(3)*ATP_concentration; 
        rates_scan(7) = scan_array(i)*ATP_concentration;

    else %Scanning something else

        rates_scan(rate_index) = scan_array(i); %Update the rate array

        %Update both ATP binding with correct ATP concentration
        rates_scan(3) = rates(3)*ATP_concentration; 
        rates_scan(7) = rates(7)*ATP_concentration;
    end

    %Calculate distribution
    [a,b] = RunCycle_OneSubunit(rates_scan, Ntrials);
    fitvals(i,:) = [a, b];


end