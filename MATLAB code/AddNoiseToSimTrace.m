function tl_out = AddNoiseToSimTrace(cfold, tl, fname, win,units)


% win = 0.1; %window size for adding noise
winsize = 1e3*win;
tl_out = tl;

for ts=1:winsize:length(tl)-winsize
    tbin = ts:ts+winsize-1;
    noise_output = GenerateNoise(cfold, win,fname, units);
    tl_out(tbin) = tl(tbin) + noise_output;

end

end