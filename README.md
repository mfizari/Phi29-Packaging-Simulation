# Phi29-Packaging-Simulation
This repository contains code for Monte Carlo simulations of DNA packaging by phi29. Code for running simulations, calculating observables, optimizing rate parameters and comparing to experimental data is written in Python (`Phi29_PackagingSimulation_Notebook`) and MATLAB (MATLAB code directory). Additionally, MATLAB code for adding experimental noise from static control tethers to simulated packaging traces is found in the MATLAB code directory. <br/> 
### Background
The purpose of these simulations is to set a threshold on the observable change in motor velocity as a function of capsid filling in [ATP] solution exchange experiments. To do this, Monte carlo simulations are used to model the mechanochemical cycle of the phi29 motor. Briefly, the motor cycles between dwell and burst phases. During a dwell, the motor's subunits sequentially exchange ADP for ATP. When all 5 subunits have bound ATP, the motor successively hydrolyzes each ATP in the burst phase, producing a 2.5 bp step for each ATP hydrolyzed except for one (the special subunit). A schematic for the dwell cycle for an invidual subunit:

![alt text](https://github.com/mfizari/Phi29-Packaging-Simulation/blob/main/Figures/s1.svg)

We also have the recently discovered "spontaneous hydrolysis" effect, where the motor can spontaneously trigger hydrolysis of <5 subunits when the dwell is long enough. This is captured by a rate constant $k_{hspont}$ (Tafoya *et al*., PNAS 2018). <br/>
### Sampling dwell times
To most efficiently simulate dwell cycles, we treat the second step in the cycle as a biased 1D random walk, where the walker can move to the $Tb$ state (tight bind) with probability $p$ or the $E$ state with probability $1-p$:

![alt text](https://github.com/mfizari/Phi29-Packaging-Simulation/blob/main/Figures/s2.svg)

The distribution of transition times to either state is exponential, with rate constants given by $k_{TB}$ for the tight binding transition and $k_{TE}$ for the ATP undocking transition. So the difference between the two transition times (whose sign tells us which state to go to), $Y=x_{T-TB} - x_{T-E}$ follows a Laplace distribution: $f(y) = e^{-yt_{TE}}/(t_{TE}+t_{TB}), y > 0$ where we have an analogous expression for $y<0$. Integreating the distribution gives: $p={k_{TE}}^{-1}/({k_{TE}}^{-1} + {k_{T-TB}}^{-1})$. Now we can easily sample a random time for the $E -> TB$ transition by drawing from the geometric distribution with parameter $p$ and finding the total time as $$\sum_{i=1}^k {t_{ET}}^{(i)} +  \sum_{i=1}^{k-1} {t_{TE}}^{(i)} + t_{TB}$$, where $k$ is the sampled number of steps and $(i)$ indicates a unique sampling from that distribution. This is implemented in the `generate_subunit_emptytotbtime` function in the Jupyter notebook. <br/>
### Rate constants
We have 7 rate constants to tune: ADP unbinding for the normal subunits and the special subunit, ATP binding for the normal subunits and the special subunit, ATP unbinding, ATP tightbinding, and spontaneous hydrolysis. We don't have to worry about hydrolysis since the rate constant for is set by measurements of burst durations (Liu *et al.*, Cell 2014). Tuning these with a grid search and minimizing the mean absolute error comparing the simulated dwell times vs [ATP] to the measured (see Jupyter notebook) sets the 7 rate constants, where we also have constraints on the relationship between the special subunit rate constants and the other subunits. Since the only rate constant that is modulated with filling is the ATP tight binding rate (Liu *et al.*, Cell 2014), we find this as a fucntion of filling by minimizing the error between measured and simulated dwell times vs filling for each measured filling level. This gives a complete set of rate constants and their filling dependence, which we compare against experimental data:

![alt text](https://github.com/mfizari/Phi29-Packaging-Simulation/blob/main/Figures/Fig1.svg)

Lastly, we can now simulate the motor vs filling for various $[ATP]$ to constrain what we would expect from measurements:

![alt text](https://github.com/mfizari/Phi29-Packaging-Simulation/blob/main/Figures/Fig2.svg)

