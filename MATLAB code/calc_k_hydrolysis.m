function kh = calc_k_hydrolysis(filling)

%This function calculates the PER SUBUNIT value of the hydrolysis rate
%Based on a linear fit of the mean burst time from 2014 Liu et al. 
%If filling < 75%, tburst = 10ms


if filling <= 75
    kh = 5/0.01;
else
    p1 = 0.002967;
    p2 = -0.2114; 
    tburst = p1*filling + p2;
    kh = 5/(tburst);

end