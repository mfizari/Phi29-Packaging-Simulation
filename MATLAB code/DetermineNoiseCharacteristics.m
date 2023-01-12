%% DetermineNoiseCharacteristics

%This script uses control tether data to determine the noise parameters

%It slides a defined window size across control tether data. In each
%window, the deviations from the mean of the window are recorded


%% Define folder directories
cfold = '';
dfold = '';
cd(dfold); files =dir('*mat'); filenames = {files.name};


%% Define parameters

window = 0.5;
slidefactor = 1; 
win = floor(1e3*window); winslide = floor(win*slidefactor);

%% Collect data

dev = [];

for i=1:length(filenames)



    %Load tether data
    cd(dfold)
    load(filenames{i}, 'l')
    cd(cfold)

    %Slide window and collect deviations
    for ts = 1:winslide:length(l)-win

        tbin = ts:ts+win;  
        lbin = l(tbin);
        lmean = mean(lbin); %mean tether length in window
        dev = [dev; lbin-lmean];

            
    end
end

dev = 1e3*dev ; %convert to bp

%% Save distribution fit from "distribution fitter" app

tagline = ['Win',num2str(window),'_Slide',num2str(slidefactor)];
fname = ['Noise_pd_',tagline,'.mat'];
save(fname, 'pd')