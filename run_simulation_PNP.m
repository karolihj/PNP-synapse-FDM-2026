% Run a simulation of the PNP model
clear all

Tstop = 10.0;      % Simulation time (ms)
dt = 0.02;         % Time step (ms)

% Define geometry (in nm)
G.Lex = 1000;
G.Lix = 600;
G.Lez = 1000;
G.Liz = 600;
G.Lm = 5;
G.Liy = 1000;
G.L_cleft = 15;
G.L_AMPA = 140;
G.Lsv = 40;
G.w_vo = 4;
G.w_ch = 4;

% Set up parameters and mesh
channels = set_up_channels(G);
G = model_parameters(G, Tstop, dt, channels);
[mesh, channels] = set_up_mesh(G, channels);

% Solve system
tic
[phi, C, states] = solve_system_PNP(G, mesh, channels);
toc

% Plot the glutamate concentration in the center of the AMPA receptor area
t = (0:G.dt:G.Tstop);
glu = reshape(C(23390,5,:), length(t), 1);
plot(t, glu)
xlabel('t (ms)')
ylabel('[Glu^-] (mM)')

