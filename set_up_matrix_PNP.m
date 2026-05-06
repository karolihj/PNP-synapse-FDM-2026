function A = set_up_matrix_PNP(G, mesh, C_prev)
%A = set_up_matrix_PNP(G, mesh, C_prev)
% Set up PDE matrix for the PNP model

% Load parameters
N = G.N;
Nx = G.Nx;
Ny = G.Ny;
elem = G.e;
kB = G.kB;
T = G.T;
dt = G.dt;
D_func = G.D_func;
eps_func = G.eps_func;


% Load grid point mesh
e_lsw = mesh.e_lsw;
e_lse = mesh.e_lse;
e_lnw = mesh.e_lnw;
e_lne = mesh.e_lne;
e_hsw = mesh.e_hsw;
e_hse = mesh.e_hse;
e_hnw = mesh.e_hnw;
e_hne = mesh.e_hne;

e_hw = mesh.e_hw;
e_he = mesh.e_he;
e_hs = mesh.e_hs;
e_hn = mesh.e_hn;
e_lw = mesh.e_lw;
e_le = mesh.e_le;
e_ls = mesh.e_ls;
e_ln = mesh.e_ln;
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;

e_w = mesh.e_w;
e_e = mesh.e_e;
e_s = mesh.e_s;
e_n = mesh.e_n;
e_h = mesh.e_h;
e_l = mesh.e_l;

e = mesh.e;

vec = zeros(N*(G.num_c+1), 1);
vec_kp = zeros(N*(G.num_c+1), 1);
vec_km = zeros(N*(G.num_c+1), 1);
vec_jp = zeros(N*(G.num_c+1), 1);
vec_jm = zeros(N*(G.num_c+1), 1);
vec_qp = zeros(N*(G.num_c+1), 1);
vec_qm = zeros(N*(G.num_c+1), 1);

vec_u = zeros(N, G.num_c);
vec_u_kp = zeros(N, G.num_c);
vec_u_km = zeros(N, G.num_c);
vec_u_jp = zeros(N, G.num_c);
vec_u_jm = zeros(N, G.num_c);
vec_u_qp = zeros(N, G.num_c);
vec_u_qm = zeros(N, G.num_c);


%% CONCENTRATION EQUATIONS
for c_idx=1:G.num_c

    zk = G.z(c_idx);
    
    %%%%%%% INNER DOMAIN %%%%%%%
    index = e;
    [x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
    dx = G.dx_values(xi);
    dy = G.dy_values(yi);
    dz = G.dz_values(zi);

    % Diffusion
    vec((c_idx-1)*N+index) = 1 + dt*(D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        D_func(x, y+dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi+1))) + ...
        D_func(x, y-dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi-1))) + ...
        D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1))));
    vec_kp((c_idx-1)*N+index+1) = -dt*D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_km((c_idx-1)*N+index-1) = -dt*D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_jp((c_idx-1)*N+index+Nx) = -dt*D_func(x, y+dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi+1)));
    vec_jm((c_idx-1)*N+index-Nx) = -dt*D_func(x, y-dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi-1)));
    vec_qp((c_idx-1)*N+index+Nx*Ny) = -dt*D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1)));
    vec_qm((c_idx-1)*N+index-Nx*Ny) = -dt*D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1)));

    % Other term
    vec_u(index, c_idx) = (dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        (dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        (dt*D_func(x, y+dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi+1))) + ...
        (dt*D_func(x, y-dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi-1))) + ...
        (dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        (dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1)));
    vec_u_kp(index+1, c_idx) = -(dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_u_km(index-1, c_idx) = -(dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_u_jp(index+Nx, c_idx) = -(dt*D_func(x, y+dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index+Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi+1))));
    vec_u_jm(index-Nx, c_idx) = -(dt*D_func(x, y-dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index-Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi-1))));
    vec_u_qp(index+Nx*Ny, c_idx) = -(dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))));
    vec_u_qm(index-Nx*Ny, c_idx) = -(dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1))));

    %%%%%%% BOUNDARY %%%%%%%

    % Neumann boundary conditions
    index = mesh.e_s_not_hold;
    [x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
    dx = G.dx_values(xi);
    dy = G.dy_values(yi);
    dz = G.dz_values(zi);
    vec((c_idx-1)*N+index) = 1 + dt*(D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        D_func(x, y+dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi+1))) + ...
        D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1))));
    vec_kp((c_idx-1)*N+index+1) = -dt*D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_km((c_idx-1)*N+index-1) = -dt*D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_jp((c_idx-1)*N+index+Nx) = -dt*D_func(x, y+dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi+1)));
    vec_qp((c_idx-1)*N+index+Nx*Ny) = -dt*D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1)));
    vec_qm((c_idx-1)*N+index-Nx*Ny) = -dt*D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1)));

    vec_u(index, c_idx) = (dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        (dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        (dt*D_func(x, y+dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi+1))) + ...
        (dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        (dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1)));
    vec_u_kp(index+1, c_idx) = -(dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_u_km(index-1, c_idx) = -(dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_u_jp(index+Nx, c_idx) = -(dt*D_func(x, y+dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index+Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi+1))));
    vec_u_qp(index+Nx*Ny, c_idx) = -(dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))));
    vec_u_qm(index-Nx*Ny, c_idx) = -(dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1))));

    index = e_n;
    [x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
    dx = G.dx_values(xi);
    dy = G.dy_values(yi);
    dz = G.dz_values(zi);
    vec((c_idx-1)*N+index) = 1 + dt*(D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        D_func(x, y-dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi-1))) + ...
        D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1))));
    vec_kp((c_idx-1)*N+index+1) = -dt*D_func(x+dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_km((c_idx-1)*N+index-1) = -dt*D_func(x-dx/2, y, z, G, c_idx)./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_jm((c_idx-1)*N+index-Nx) = -dt*D_func(x, y-dy/2, z, G, c_idx)./(dy.*0.5.*(dy+G.dy_values(yi-1)));
    vec_qp((c_idx-1)*N+index+Nx*Ny) = -dt*D_func(x, y, z+dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi+1)));
    vec_qm((c_idx-1)*N+index-Nx*Ny) = -dt*D_func(x, y, z-dz/2, G, c_idx)./(dz.*0.5.*(dz+G.dz_values(zi-1)));

    vec_u(index, c_idx) = (dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1))) + ...
        (dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1))) + ...
        (dt*D_func(x, y-dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi-1))) + ...
        (dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))) + ...
        (dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1)));
    vec_u_kp(index+1, c_idx) = -(dt*D_func(x+dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index+1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi+1)));
    vec_u_km(index-1, c_idx) = -(dt*D_func(x-dx/2, y, z, G, c_idx)*zk*elem/(2*kB*T)).*(C_prev(index-1,c_idx) + C_prev(index,c_idx))./(dx.*0.5.*(dx+G.dx_values(xi-1)));
    vec_u_jm(index-Nx, c_idx) = -(dt*D_func(x, y-dy/2, z, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index-Nx,c_idx) + C_prev(index,c_idx))./(dy.*0.5.*(dy+G.dy_values(yi-1))));
    vec_u_qp(index+Nx*Ny, c_idx) = -(dt*D_func(x, y, z+dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index+Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi+1))));
    vec_u_qm(index-Nx*Ny, c_idx) = -(dt*D_func(x, y, z-dz/2, G, c_idx)*zk*elem/(2*kB*T)).*((C_prev(index-Nx*Ny,c_idx) + C_prev(index,c_idx))./(dz.*0.5.*(dz+G.dz_values(zi-1))));

    % Dirichlet boundary conditions
    index = [e_h, e_l, e_w, e_e, e_nw, e_sw, e_se, e_ne, e_ln, e_hn, ...
        e_ls, e_hs, e_hw, e_lw, e_he, e_le, e_lsw, e_lse, e_lnw, ...
        e_lne, e_hsw, e_hse, e_hnw, e_hne, mesh.hold'];
    vec((c_idx-1)*N+index) = 1;
        
end


%% PHI EQUATION

%%%%%%% INNER DOMAIN %%%%%%%
index = mesh.e;
[x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
dx = G.dx_values(xi);
dy = G.dy_values(yi);
dz = G.dz_values(zi);
vec(G.num_c*N+index) = -(eps_func(x+dx/2, y, z, G)./(dx.*0.5.*(dx+G.dx_values(xi+1))) ...
    + eps_func(x-dx/2, y, z, G)./(dx*0.5.*(dx+G.dx_values(xi-1))) ...
    + eps_func(x, y+dy/2, z, G)./(dy.*0.5.*(dy+G.dy_values(yi+1))) ...
    + eps_func(x, y-dy/2, z, G)./(dy.*0.5.*(dy+G.dy_values(yi-1))) ...
    + eps_func(x, y, z+dz/2, G)./(dz.*0.5.*(dz+G.dz_values(zi+1))) ...
    + eps_func(x, y, z-dz/2, G)./(dz.*0.5.*(dz+G.dz_values(zi-1))));
vec_kp(G.num_c*N+index+1) = eps_func(x+dx/2, y, z, G)./(dx.*0.5.*(dx+G.dx_values(xi+1)));
vec_km(G.num_c*N+index-1) = eps_func(x-dx/2, y, z, G)./(dx.*0.5.*(dx+G.dx_values(xi-1)));
vec_jp(G.num_c*N+index+Nx) = eps_func(x, y+dy/2, z, G)./(dy.*0.5.*(dy+G.dy_values(yi+1)));
vec_jm(G.num_c*N+index-Nx) = eps_func(x, y-dy/2, z, G)./(dy.*0.5.*(dy+G.dy_values(yi-1)));
vec_qp(G.num_c*N+index+Nx*Ny) = eps_func(x, y, z+dz/2, G)./(dz.*0.5.*(dz+G.dz_values(zi+1)));
vec_qm(G.num_c*N+index-Nx*Ny) = eps_func(x, y, z-dz/2, G)./(dz.*0.5.*(dz+G.dz_values(zi-1)));

%%%%%%% BOUNDARY %%%%%%%

% Dirichlet boundary condition
index = [e_h, e_l, e_w, e_e, e_nw, e_sw, e_se, e_ne, e_ln, e_hn, e_ls, e_hs, ...
        e_hw, e_lw, e_he, e_le, e_lsw, e_lse, e_lnw, e_lne, e_hsw, e_hse, ...
        e_hnw, e_hne, mesh.hold'];
vec(G.num_c*N+index) = 1;

% Neumann boundary condition
index = e_n;
[x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
dx = G.dx_values(xi);
dy = G.dy_values(yi);
dz = G.dz_values(zi);
vec(G.num_c*N+index) = -eps_func(x, y-dy/2, z, G)./dy;
vec_jm(G.num_c*N+index-Nx) = eps_func(x, y-dy/2, z, G)./dy;

index = mesh.e_s_not_hold;
[x, y, z, xi, yi, zi] = x_y_z_from_idx(index, G);
dx = G.dx_values(xi);
dy = G.dy_values(yi);
dz = G.dz_values(zi);
vec(G.num_c*N+index) = -eps_func(x, y+dy/2, z, G)./dy;
vec_jp(G.num_c*N+index+Nx) = eps_func(x, y+dy/2, z, G)./dy;



%%%%%%% SET UP THE MATRIX %%%%%%%
A = spdiags(vec, 0, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_kp, 1, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_km, -1, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_jp, Nx, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_jm, -Nx, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_qp, Nx*Ny, N*(G.num_c+1), N*(G.num_c+1));
A = A + spdiags(vec_qm, -Nx*Ny, N*(G.num_c+1), N*(G.num_c+1));

for c_idx=1:G.num_c
    A_tmp = spdiags(vec_u(:,c_idx), 0, N, N);
    A_tmp = A_tmp + spdiags(vec_u_kp(:,c_idx), 1, N, N);
    A_tmp = A_tmp + spdiags(vec_u_km(:,c_idx), -1, N, N);
    A_tmp = A_tmp + spdiags(vec_u_jp(:,c_idx),Nx, N, N);
    A_tmp = A_tmp + spdiags(vec_u_jm(:,c_idx), -Nx, N, N);
    A_tmp = A_tmp + spdiags(vec_u_qp(:,c_idx), Nx*Ny, N, N);
    A_tmp = A_tmp + spdiags(vec_u_qm(:,c_idx), -Nx*Ny, N, N);
    A((c_idx-1)*N+(1:N), G.num_c*N+(1:N)) = A((c_idx-1)*N+(1:N), G.num_c*N+(1:N)) + A_tmp;
end

% CONCENTRATIONS IN PHI
vec = zeros(G.N, 1);
vec(mesh.e) = 1;
for c_idx = 1:G.num_c
    A_tmp = spdiags(G.F*G.z(c_idx)*vec, 0, N, N);
    A(G.num_c*N+(1:N), (c_idx-1)*N+(1:N)) = A_tmp;
end


end

