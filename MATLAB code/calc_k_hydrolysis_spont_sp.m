function k_hydrolysis_spont_sp = calc_k_hydrolysis_spont_sp(filling)

   xdata = [0,100]; 
   ydata = [0.4*4, 0.4*4];

   p = polyfit(xdata,ydata,1);


    k_hydrolysis_spont_sp = p(1)*filling + p(2);



end