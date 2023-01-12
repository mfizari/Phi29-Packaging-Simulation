function k_ATP_tightbind = calc_k_ATP_tightbind(filling)

%This function calculates a k_ATP_tightbind rate given a filling level

fillplot = [0;20;32.5;55;65;75;87.5;92.5;97.5]';
k_ATP_tightbind_fitdata = 1e3*[1.333,1.333, 1.333,0.088, 0.0451, 0.0308, 0.0125, 0.0088, 0.0045]';

x2 = fillplot(fillplot >= 32.5);
y2 = k_ATP_tightbind_fitdata(fillplot >= 32.5);

snerbfactor = 1;

generate = fit(x2',snerbfactor*y2,'Power1');


if filling <=32.5
    k_ATP_tightbind = k_ATP_tightbind_fitdata(1);
else
    k_ATP_tightbind = generate(filling);
end





