function channels = set_up_channels(G)
%channels = set_up_channels(G)
% Note: only AMPA receptors (channels{13}) are open in the simulations

G.Lx = 2*G.Lex + 2*G.Lm + G.Lix;
G.Ly = 2*G.Liy + 2*G.Lm + G.L_cleft;
G.Lz = 2*G.Lez + 2*G.Lm + G.Liz;

% Potassium channel presynaptic cell (left)
channels{1}.lx = G.Lm;         % nm 
channels{1}.ly = G.w_ch*5;     % nm 
channels{1}.lz = G.w_ch*5;     % nm  
channels{1}.coordinate = 'x';
channels{1}.direction = 'ei';
channels{1}.A = channels{1}.ly*channels{1}.lz;
channels{1}.channel_start_x = G.Lex;
channels{1}.channel_start_y = G.Ly - (G.Liy/8 + channels{1}.ly/2);
channels{1}.channel_start_z = G.Lz/2-channels{1}.lz/2;
channels{1}.g = 5;             % pS 
channels{1}.nc = 25;
channels{1}.flux = @ion_channel_flux;
channels{1}.ion = 2;
channels{1}.open = @(t) 0; 
channels{1}.num_states = 0;
channels{1}.rhs_states = @AMPA_receptor_kinetics;

% Potassium channel presynaptic cell (right)
channels{2}.lx = G.Lm;         % nm 
channels{2}.ly = G.w_ch*5;     % nm 
channels{2}.lz = G.w_ch*5;     % nm
channels{2}.coordinate = 'x';
channels{2}.direction = 'ie';
channels{2}.A = channels{2}.ly*channels{2}.lz;
channels{2}.channel_start_x = G.Lex + G.Lm + G.Lix;
channels{2}.channel_start_y = G.Ly - (G.Liy/8 + channels{2}.ly/2);
channels{2}.channel_start_z = G.Lz/2-channels{2}.lz/2;
channels{2}.g = 5;             % pS 
channels{2}.nc = 25;
channels{2}.flux = @ion_channel_flux;
channels{2}.ion = 2;
channels{2}.open = @(t) 0; 
channels{2}.num_states = 0;
channels{2}.rhs_states = @AMPA_receptor_kinetics;


% Potassium channel postsynaptic cell (left)
channels{3}.lx = G.Lm;         % nm 
channels{3}.ly = G.w_ch*5;     % nm 
channels{3}.lz = G.w_ch*5;     % nm   
channels{3}.coordinate = 'x';
channels{3}.direction = 'ei';
channels{3}.A = channels{3}.ly*channels{3}.lz;
channels{3}.channel_start_x = G.Lex;
channels{3}.channel_start_y = G.Liy/8 - channels{3}.ly/2;
channels{3}.channel_start_z = G.Lz/2-channels{3}.lz/2;
channels{3}.g = 5;             % pS 
channels{3}.nc = 25;
channels{3}.flux = @ion_channel_flux;
channels{3}.ion = 2;
channels{3}.open = @(t) 0; 
channels{3}.num_states = 0;
channels{3}.rhs_states = @AMPA_receptor_kinetics;

% Potassium channel postsynaptic cell (right)
channels{4}.lx = G.Lm;         % nm 
channels{4}.ly = G.w_ch*5;     % nm 
channels{4}.lz = G.w_ch*5;     % nm  
channels{4}.coordinate = 'x';
channels{4}.direction = 'ie';
channels{4}.A = channels{4}.ly*channels{4}.lz;
channels{4}.channel_start_x = G.Lex + G.Lm + G.Lix;
channels{4}.channel_start_y = G.Liy/8 - channels{4}.ly/2;
channels{4}.channel_start_z = G.Lz/2-channels{4}.lz/2;
channels{4}.g = 5;             % pS 
channels{4}.nc = 25;
channels{4}.flux = @ion_channel_flux;
channels{4}.ion = 2;
channels{4}.open = @(t) 0; 
channels{4}.num_states = 0;
channels{4}.rhs_states = @AMPA_receptor_kinetics;

% Sodium channel presynaptic cell (left)
channels{5}.lx = G.Lm;         % nm 
channels{5}.ly = G.w_ch*5;     % nm 
channels{5}.lz = G.w_ch*5;     % nm  
channels{5}.coordinate = 'x';
channels{5}.direction = 'ei';
channels{5}.A = channels{5}.ly*channels{5}.lz;
channels{5}.channel_start_x = G.Lex;
channels{5}.channel_start_y = G.Ly - (G.Liy/4 + channels{5}.ly/2);
channels{5}.channel_start_z = G.Lz/2-channels{5}.lz/2;
channels{5}.g = 20;             % pS 
channels{5}.nc = 25;
channels{5}.flux = @ion_channel_flux;
channels{5}.ion = 1;
channels{5}.open = @(t) 0; 
channels{5}.num_states = 0;
channels{5}.rhs_states = @AMPA_receptor_kinetics;

% Sodium channel presynaptic cell (right)
channels{6}.lx = G.Lm;         % nm 
channels{6}.ly = G.w_ch*5;     % nm 
channels{6}.lz = G.w_ch*5;     % nm   
channels{6}.coordinate = 'x';
channels{6}.direction = 'ie';
channels{6}.A = channels{6}.ly*channels{6}.lz;
channels{6}.channel_start_x = G.Lex + G.Lm + G.Lix;
channels{6}.channel_start_y = G.Ly - (G.Liy/4 + channels{6}.ly/2);
channels{6}.channel_start_z = G.Lz/2-channels{6}.lz/2;
channels{6}.g = 20;            % pS 
channels{6}.nc = 25;
channels{6}.flux = @ion_channel_flux;
channels{6}.ion = 1;
channels{6}.open = @(t) 0; 
channels{6}.num_states = 0;
channels{6}.rhs_states = @AMPA_receptor_kinetics;


% Sodium channel postsynaptic cell (left)
channels{7}.lx = G.Lm;         % nm 
channels{7}.ly = G.w_ch*5;     % nm 
channels{7}.lz = G.w_ch*5;     % nm  
channels{7}.coordinate = 'x';
channels{7}.direction = 'ei';
channels{7}.A = channels{7}.ly*channels{7}.lz;
channels{7}.channel_start_x = G.Lex;
channels{7}.channel_start_y = G.Liy/4 - channels{7}.ly/2;
channels{7}.channel_start_z = G.Lz/2-channels{7}.lz/2;
channels{7}.g = 20;             % pS 
channels{7}.nc = 25;
channels{7}.flux = @ion_channel_flux;
channels{7}.ion = 1;
channels{7}.open = @(t) 0; 
channels{7}.num_states = 0;
channels{7}.rhs_states = @AMPA_receptor_kinetics;


% Sodium channel postsynaptic cell (right)
channels{8}.lx = G.Lm;         % nm 
channels{8}.ly = G.w_ch*5;     % nm 
channels{8}.lz = G.w_ch*5;     % nm 
channels{8}.coordinate = 'x';
channels{8}.direction = 'ie';
channels{8}.A = channels{8}.ly*channels{8}.lz;
channels{8}.channel_start_x = G.Lex + G.Lm + G.Lix;
channels{8}.channel_start_y = G.Liy/8 - channels{8}.ly/2;
channels{8}.channel_start_z = G.Lz/2-channels{8}.lz/2;
channels{8}.g = 20;            % pS 
channels{8}.nc = 25;
channels{8}.flux = @ion_channel_flux;
channels{8}.ion = 1;
channels{8}.open = @(t) 0; 
channels{8}.num_states = 0;
channels{8}.rhs_states = @AMPA_receptor_kinetics;


% Calcium channel presynaptic cell (left)
channels{9}.lx = G.w_ch;     % nm 
channels{9}.ly = G.Lm;       % nm 
channels{9}.lz = G.w_ch;     % nm 
channels{9}.coordinate = 'y';
channels{9}.direction = 'ei';
channels{9}.A = channels{9}.lx*channels{9}.lz;
channels{9}.channel_start_x = G.Lex + G.Lm + (G.Lix-G.Lsv)/2 - G.w_ch - channels{9}.lx; 
channels{9}.channel_start_y = G.Liy + G.Lm + G.L_cleft;
channels{9}.channel_start_z = G.Lz/2-channels{9}.lz/2;
channels{9}.g = 2;            % pS 
channels{9}.nc = 1;
channels{9}.flux = @ion_channel_flux;
channels{9}.ion = 3;
channels{9}.open = @(t) 0; 
channels{9}.num_states = 0;
channels{9}.rhs_states = @AMPA_receptor_kinetics;


% Calcium channel presynaptic cell (right)
channels{10}.lx = G.w_ch;     % nm 
channels{10}.ly = G.Lm;       % nm 
channels{10}.lz = G.w_ch;     % nm 
channels{10}.coordinate = 'y';
channels{10}.direction = 'ei';
channels{10}.A = channels{10}.lx*channels{10}.lz;
channels{10}.channel_start_x = G.Lex + G.Lm + (G.Lix+G.Lsv)/2 + G.w_ch; 
channels{10}.channel_start_y = G.Liy + G.Lm + G.L_cleft;
channels{10}.channel_start_z = G.Lz/2-channels{10}.lz/2;
channels{10}.g = 2;           % pS 
channels{10}.nc = 1;
channels{10}.flux = @ion_channel_flux;
channels{10}.ion = 3;
channels{10}.open = @(t) 0; 
channels{10}.num_states = 0;
channels{10}.rhs_states = @AMPA_receptor_kinetics;

% Calcium channel presynaptic cell (lower)
channels{11}.lx = G.w_ch;     % nm 
channels{11}.ly = G.Lm;       % nm 
channels{11}.lz = G.w_ch;     % nm  
channels{11}.coordinate = 'y';
channels{11}.direction = 'ei';
channels{11}.A = channels{11}.lx*channels{11}.lz;
channels{11}.channel_start_x = G.Lx/2-channels{11}.lx/2;
channels{11}.channel_start_y = G.Liy + G.Lm + G.L_cleft;
channels{11}.channel_start_z = G.Lez + G.Lm + (G.Liz-G.Lsv)/2 - G.w_ch - channels{11}.lz; 
channels{11}.g = 2;           % pS 
channels{11}.nc = 1;
channels{11}.flux = @ion_channel_flux;
channels{11}.ion = 3;
channels{11}.open = @(t) 0; 
channels{11}.num_states = 0;
channels{11}.rhs_states = @AMPA_receptor_kinetics;


% Calcium channel presynaptic cell (upper)
channels{12}.lx = G.w_ch;     % nm 
channels{12}.ly = G.Lm;       % nm 
channels{12}.lz = G.w_ch;     % nm   
channels{12}.coordinate = 'y';
channels{12}.direction = 'ei';
channels{12}.A = channels{12}.lx*channels{12}.lz;
channels{12}.channel_start_x = G.Lx/2-channels{12}.lx/2;
channels{12}.channel_start_y = G.Liy + G.Lm + G.L_cleft;
channels{12}.channel_start_z = G.Lez + G.Lm + (G.Liz+G.Lsv)/2 + G.w_ch; 
channels{12}.g = 2;           % pS 
channels{12}.nc = 1;
channels{12}.flux = @ion_channel_flux;
channels{12}.ion = 3;
channels{12}.open = @(t) 0; 
channels{12}.num_states = 0;
channels{12}.rhs_states = @AMPA_receptor_kinetics;


% AMPA receptors postsynaptic cell
channels{13}.lx = G.L_AMPA;     % nm 
channels{13}.ly = G.Lm;         % nm 
channels{13}.lz = G.L_AMPA;     % nm  
channels{13}.coordinate = 'y';
channels{13}.direction = 'ie';
channels{13}.A = channels{13}.lx*channels{13}.lz;
channels{13}.channel_start_x = G.Lex + G.Lm + (G.Lix-G.L_AMPA)/2;
channels{13}.channel_start_y = G.Liy;
channels{13}.channel_start_z = G.Lez + G.Lm + (G.Liz-G.L_AMPA)/2;
channels{13}.g = [25, 15];      % pS
channels{13}.nc = 200;
channels{13}.flux = @ion_channel_flux;
channels{13}.ion = [1; 2];
channels{13}.open = @(t) 0; 
channels{13}.num_states = 6;
channels{13}.rhs_states = @AMPA_receptor_kinetics;


end
