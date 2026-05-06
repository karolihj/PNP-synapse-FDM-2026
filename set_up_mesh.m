function [mesh, channels] = set_up_mesh(G, channels)
%[mesh, channels] = set_up_mesh(G, channels) 
% Set up vectors containing the indices of the different types of grid nodes 
% and prepare channels

% Read geometry
Nx = G.Nx;
Ny = G.Ny;
Nz = G.Nz;
N = G.N;

% Set up indices for outer boundary
% Set up indices for outer boundary (corners)
mesh.e_lsw = 1;
mesh.e_lse = Nx;
mesh.e_lnw = (Ny-1)*Nx + 1;
mesh.e_lne = Nx*Ny;
mesh.e_hsw = (Nz-1)*Ny*Nx+1;
mesh.e_hse = (Nz-1)*Ny*Nx+Nx;
mesh.e_hnw = ((Nz-1)*Ny+Ny-1)*Nx + 1;
mesh.e_hne = Nx*Ny*Nz;

% Set up indices for outer boundary (lines)
mesh.e_lw = ((2:Ny-1)-1)*Nx + 1;
mesh.e_le = ((2:Ny-1)-1)*Nx + Nx;
mesh.e_ls = 2:Nx-1;
mesh.e_ln = (Ny-1)*Nx + (2:Nx-1);
mesh.e_hw = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + 1;
mesh.e_he = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + Nx;
mesh.e_hs = (Nz-1)*Ny*Nx+(2:Nx-1);
mesh.e_hn = ((Nz-1)*Ny+Ny-1)*Nx + (2:Nx-1);
mesh.e_sw = (1:Nz-2)*Ny*Nx+1;
mesh.e_se = (1:Nz-2)*Ny*Nx+Nx;
mesh.e_nw = ((1:Nz-2)*Ny+Ny-1)*Nx+1;
mesh.e_ne = ((1:Nz-2)*Ny+Ny-1)*Nx+Nx;

% Set up indices for outer boundary (sides)
[x, y] = meshgrid(2:Nx-1, 2:Ny-1);
mesh.e_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
mesh.e_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, Nz*ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
[x, z] = meshgrid(2:Nx-1, 2:Nz-1);
mesh.e_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
mesh.e_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, Ny*ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
[y, z] = meshgrid(2:Ny-1, 2:Nz-1);
mesh.e_w = sort(reshape(sub2ind([Nx,Ny,Nz], ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';
mesh.e_e =  sort(reshape(sub2ind([Nx,Ny,Nz], Nx*ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';


% Set up channel indices
% Channel indices
[x, y, z] = x_y_z_from_idx(1:G.N, G);
for i=1:length(channels)
    if strcmp(channels{i}.coordinate, 'x') && strcmp(channels{i}.direction, 'ei')
        channels{i}.m_e = find((x > channels{i}.channel_start_x-G.dx_min/2-G.dx_min/2).*(x < channels{i}.channel_start_x-G.dx_min/2+G.dx_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
        channels{i}.m_i = find((x > channels{i}.channel_start_x+G.Lm+G.dx_min/2-G.dx_min/4).*(x < channels{i}.channel_start_x+G.Lm+G.dx_min/2+G.dx_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
    elseif strcmp(channels{i}.coordinate, 'x') && strcmp(channels{i}.direction, 'ie')
        channels{i}.m_i = find((x > channels{i}.channel_start_x-G.dx_min/2-G.dx_min/2).*(x < channels{i}.channel_start_x-G.dx_min/2+G.dx_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
        channels{i}.m_e = find((x > channels{i}.channel_start_x+G.Lm+G.dx_min/2-G.dx_min/4).*(x < channels{i}.channel_start_x+G.Lm+G.dx_min/2+G.dx_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
    elseif strcmp(channels{i}.coordinate, 'y') && strcmp(channels{i}.direction, 'ei')
        channels{i}.m_e = find((y > channels{i}.channel_start_y-G.dy_min/2-G.dy_min/2).*(y < channels{i}.channel_start_y-G.dy_min/2+G.dy_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx));
        channels{i}.m_i = find((y > channels{i}.channel_start_y+G.Lm+G.dy_min/2-G.dy_min/4).*(y < channels{i}.channel_start_y+G.Lm+G.dy_min/2+G.dy_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx));
    elseif strcmp(channels{i}.coordinate, 'y') && strcmp(channels{i}.direction, 'ie')
        channels{i}.m_i = find((y > channels{i}.channel_start_y-G.dy_min/2-G.dy_min/2).*(y < channels{i}.channel_start_y-G.dy_min/2+G.dy_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx));
        channels{i}.m_e = find((y > channels{i}.channel_start_y+G.Lm+G.dy_min/2-G.dy_min/4).*(y < channels{i}.channel_start_y+G.Lm+G.dy_min/2+G.dy_min/4).*...
            (z >= channels{i}.channel_start_z).*(z <= channels{i}.channel_start_z+channels{i}.lz).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx));
    elseif strcmp(channels{i}.coordinate, 'z') && strcmp(channels{i}.direction, 'ei')
        channels{i}.m_e = find((z > channels{i}.channel_start_z-G.dz_min/2-G.dz_min/2).*(z < channels{i}.channel_start_z-G.dz_min/2+G.dz_min/4).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
        channels{i}.m_i = find((z > channels{i}.channel_start_z+G.Lm+G.dz_min/2-G.dz_min/4).*(z < channels{i}.channel_start_z+G.Lm+G.dz_min/2+G.dz_min/4).*...
            (x >= channels{i}.channel_start_x).*(z <= channels{i}.channel_start_x+channels{i}.lx).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
    elseif strcmp(channels{i}.coordinate, 'z') && strcmp(channels{i}.direction, 'ie')
        channels{i}.m_i = find((z > channels{i}.channel_start_z-G.dz_min/2-G.dz_min/2).*(z < channels{i}.channel_start_z-G.dz_min/2+G.dz_min/4).*...
            (x >= channels{i}.channel_start_x).*(x <= channels{i}.channel_start_x+channels{i}.lx).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
        channels{i}.m_e = find((z > channels{i}.channel_start_z+G.Lm+G.dz_min/2-G.dz_min/4).*(z < channels{i}.channel_start_z+G.Lm+G.dz_min/2+G.dz_min/4).*...
            (x >= channels{i}.channel_start_x).*(z <= channels{i}.channel_start_x+channels{i}.lx).*...
            (y >= channels{i}.channel_start_y).*(y <= channels{i}.channel_start_y+channels{i}.ly));
    end

    % State variables
    channels{i}.states = zeros(length(channels{i}.m_i), channels{i}.num_states);

end

% Set up initial conditions for the state variables
channels{13}.states(:,1) = 1; % Correct initial conditions

% Set up indices for holding potential
mesh.hold  = find((y < G.dy_min).*(z >= G.Lz/2-G.L_AMPA/2).*(z <= G.Lz/2+G.L_AMPA/2).*...
    (x >= G.Lx/2-G.L_AMPA/2).*(x <= G.Lx/2+G.L_AMPA/2));
all = zeros(G.N,1);
all(mesh.e_s) = 1;
all(mesh.hold) = 0;
mesh.e_s_not_hold = find(all);


% Set up indices for the remaining voxels
indices = ones(1, N);
indices([mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse, ...
    mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, ...
    mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, ...
    mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
    mesh.e_l]) = 0;
indices(mesh.hold) = 0;
mesh.e = round(find(indices));

end
