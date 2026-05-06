function G = model_parameters(G, Tstop, dt, channels)
%G = model_parameters(G, Tstop, dt, channels)

G.Lx = 2*G.Lex + 2*G.Lm + G.Lix;
G.Ly = 2*G.Liy + 2*G.Lm + G.L_cleft;
G.Lz = 2*G.Lez + 2*G.Lm + G.Liz;

G.dt = dt;
G.Nt = round(Tstop/dt);
G.Tstop = Tstop;
G.dt_ode = 0.001;
G.nt = round(G.dt/G.dt_ode);

% Set up mesh
G.dx_min = 5;     % Near membrane
G.dx_ch = 4;      % Near channels
G.dx_AMPA = 14;   % AMPA receptor area
G.dx_max = 10000;
G.dy_min = 1;     % Near membrane
G.dy_ch = 20;     % Near channels
G.dy_max = 10000;
G.dz_min = 5;     % Near membrane
G.dz_ch = 4;      % Near channels
G.dz_AMPA = 14;   % AMPA receptor area
G.dz_max = 10000;
G.num_boundary_layer_dx = 1;
G.increase_factor_x = 8;
G.num_per_dx_increase = 1;
G.num_boundary_layer_dy = 2;
G.increase_factor_y = 4;
G.num_per_dy_increase = 1;
G.num_boundary_layer_dz = 1;
G.increase_factor_z = 8;
G.num_per_dz_increase = 1;

[G.dx_values, G.Nx] = set_up_dx(G);
[G.dy_values, G.Ny] = set_up_dy(G, channels);
[G.dz_values, G.Nz] = set_up_dz(G);
G.N = G.Nx*G.Ny*G.Nz;

% Physical constants
G.num_c = 5; % (Na, K, Ca, Cl, glut)
G.z = [1, 1, 2, -1, -1];
G.D = [1.33e6; 1.96e6; 0.71e6; 2.03e6; 0.86e6]; % nm^2/ms
G.kappa_cleft = [2.56; 2.56; 2.56; 2.56; 2.56];
G.D_func = @diffusion_coefficient;
G.c0e = [100; 4; 1.4; 106.8; 0]; % mM
G.c0i = [12; 125; 0.0001; 137.0002; 0]; % mM
G.c0v = [145; 5; 0; 0; 150]; % mM 
G.rho0 = zeros(G.N, 1);
G.initial_conditions = @initial_conditions;

G.eps_r = 80;
G.eps_m = 2;
G.eps_c = 2;
G.eps_0 = 8854; %fF/m 
G.eps = G.eps_r*G.eps_0;
G.eps_mem = G.eps_m*G.eps_0;
G.eps_channel = G.eps_r*G.eps_0;
G.eps_func = @permittivity;

G.e = 1.60217662e-19; % C
G.kB = 1.38064852e-20; % mJ/K
G.T = 310; % K (300)
G.NA = 6.02214076e20; % mmol^{-1}
G.F = 96485.3365; % C/mol

% Holding potential
G.V_hold = -70;

% Open vesicle
G.open_vesicle = 1;

end

