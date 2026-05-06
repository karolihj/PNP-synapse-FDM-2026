function [phi, C, states] = solve_system_PNP(G, mesh, channels)
%[phi, C, states] = solve_system_PNP(G, mesh, channels)
% Run a simulation of the PNP model

% Set up initial conditions
[x, y, z] = x_y_z_from_idx(1:G.N, G);
c = zeros(G.N, G.num_c);
for i=1:G.num_c
    c(:, i) = G.initial_conditions(x, y, z, G, i);
end


% Set up matrices for saving the solution
G.num_save = G.Nt+1;
phi = zeros(G.N, G.num_save);
U = zeros(G.N, 1);
C = zeros(G.N, G.num_c, G.num_save);
C(:,:,1) = c;
t = 0;

states = zeros(size(channels{13}.states, 1), size(channels{13}.states, 2), G.num_save);
states(:,:,1) = channels{13}.states;


% Run simulation
fprintf('Starting simulation...\n')
print_step = G.dt;
n_print = max(round(print_step/G.dt), 1);
t1 = tic;
for n = 1:G.Nt
    
    t = t + G.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 1: Solve AMPA receptor kinetics model   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:length(channels)
        if channels{i}.num_states > 0
            for k=1:G.nt
                channels{i}.states = channels{i}.states + G.dt_ode*channels{i}.rhs_states(channels{i}.states, c(channels{i}.m_e,:));
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 2: Solve PNP model               %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bc = c;

    % Set up matrix
    A = set_up_matrix_PNP(G, mesh, c);

    % Update channel fluxes
    for i=1:length(channels)
        if strcmp(channels{i}.coordinate, 'x')
            dx = G.dx_min;
        elseif strcmp(channels{i}.coordinate, 'y')
            dx = G.dy_min;
        else
            dx = G.dz_min;
        end
        for k=1:length(channels{i}.ion)
            ion = channels{i}.ion(k);
            [J, J_rhs, J_factor] = channels{i}.flux(c, U, G, channels{i}, k, t);
            if channels{i}.num_states > 0
                ch_open = 1 - sum(channels{i}.states,2);
            else
                ch_open = channels{i}.open(t)*ones(size(channels{i}.m_e));
            end
            for j=1:length(channels{i}.m_e)
                if ch_open(j) > 0
                    bc(channels{i}.m_e(j), ion) = bc(channels{i}.m_e(j), ion) + G.dt*J_rhs(j)/dx;
                    bc(channels{i}.m_i(j), ion) = bc(channels{i}.m_i(j), ion) - G.dt*J_rhs(j)/dx;

                    A((ion-1)*G.N+channels{i}.m_e(j), G.num_c*G.N+channels{i}.m_e(j)) = A((ion-1)*G.N+channels{i}.m_e(j), G.num_c*G.N+channels{i}.m_e(j)) + G.dt*J_factor(j)/dx;
                    A((ion-1)*G.N+channels{i}.m_e(j), G.num_c*G.N+channels{i}.m_i(j)) = A((ion-1)*G.N+channels{i}.m_e(j), G.num_c*G.N+channels{i}.m_i(j)) - G.dt*J_factor(j)/dx;
                    A((ion-1)*G.N+channels{i}.m_i(j), G.num_c*G.N+channels{i}.m_e(j)) = A((ion-1)*G.N+channels{i}.m_i(j), G.num_c*G.N+channels{i}.m_e(j)) - G.dt*J_factor(j)/dx;
                    A((ion-1)*G.N+channels{i}.m_i(j), G.num_c*G.N+channels{i}.m_i(j)) = A((ion-1)*G.N+channels{i}.m_i(j), G.num_c*G.N+channels{i}.m_i(j)) + G.dt*J_factor(j)/dx;

                    A((ion-1)*G.N+channels{i}.m_e(j), (ion-1)*G.N+channels{i}.m_e(j)) = A((ion-1)*G.N+channels{i}.m_e(j), (ion-1)*G.N+channels{i}.m_e(j)) + G.dt*(J_factor(j)*(G.kB*G.T/(G.z(ion)*G.e))./c(channels{i}.m_e(j), ion))/dx;
                    A((ion-1)*G.N+channels{i}.m_e(j), (ion-1)*G.N+channels{i}.m_i(j)) = A((ion-1)*G.N+channels{i}.m_e(j), (ion-1)*G.N+channels{i}.m_i(j)) - G.dt*(J_factor(j)*(G.kB*G.T/(G.z(ion)*G.e))./c(channels{i}.m_i(j), ion))/dx;
                    A((ion-1)*G.N+channels{i}.m_i(j), (ion-1)*G.N+channels{i}.m_e(j)) = A((ion-1)*G.N+channels{i}.m_i(j), (ion-1)*G.N+channels{i}.m_e(j)) - G.dt*(J_factor(j)*(G.kB*G.T/(G.z(ion)*G.e))./c(channels{i}.m_e(j), ion))/dx;
                    A((ion-1)*G.N+channels{i}.m_i(j), (ion-1)*G.N+channels{i}.m_i(j)) = A((ion-1)*G.N+channels{i}.m_i(j), (ion-1)*G.N+channels{i}.m_i(j)) + G.dt*(J_factor(j)*(G.kB*G.T/(G.z(ion)*G.e))./c(channels{i}.m_i(j), ion))/dx;
                end
            end
        end
    end


    % Update rhs for phi equation
    b_phi = zeros(G.N, 1);
    b_phi(mesh.e) = -G.rho0(mesh.e);
    b_phi(mesh.hold) = G.V_hold;

    % Solve system
    b_tot = [reshape(bc, G.num_c*G.N, 1); b_phi];
    X = A\b_tot;

    % Extract solutions
    c = reshape(X(1:G.N*G.num_c), G.N, G.num_c);
    U = X(G.N*G.num_c+1:end);

    % Save the solution
    phi(:,n+1) = U;
    C(:,:,n+1) = c;
    states(:,:,n+1) = channels{13}.states;
    
    % Estimate remaining simulation time
    if rem(n, n_print) == 0
        % Print current point in time
        fprintf('t = %.3g ms. ', t);

        % Print estimated simulation time
        t2 = toc(t1);                % Time usage for n_print time steps
        t_rem = t2*(G.Nt-n)/n_print; % Estimated remaining simulation time
        fprintf('EstRemSimTime: ');
        if t_rem > 86400*2
            fprintf('%.1f days \n', t_rem/86400);
        elseif t_rem > 3600
            fprintf('%.1f h \n', t_rem/3600);
        elseif t_rem > 60
            fprintf('%.1f min \n', t_rem/60);
        else
            fprintf('%.1f sec \n', t_rem);
        end
        t1 = tic;
    end
end

end

