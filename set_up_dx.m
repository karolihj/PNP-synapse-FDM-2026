function [dx_values, Nx] = set_up_dx(G)
% [dx_values, Nx] = set_up_dx(G)
% Set up adaptive mesh in the x-direction

increase_factor = G.increase_factor_x;

%% Extracellular space
dx_valuesE = G.dx_min*ones(G.num_boundary_layer_dx, 1);

x = G.num_boundary_layer_dx*G.dx_min + G.dx_min/2;
counter = 0;
if G.dx_max >= G.dx_min*increase_factor
    current_dx = G.dx_min*increase_factor;
else
    current_dx = G.dx_min;
end
while x < G.Lex -G.dx_min
    dx_valuesE = [dx_valuesE; current_dx];
    x = x + current_dx;
    counter = counter + 1;
    if counter == G.num_per_dx_increase
        if current_dx*increase_factor <= G.dx_max
            current_dx = current_dx*increase_factor;
        end
        counter = 0;
    end
end

dx_valuesE(end) = G.Lex -(sum(dx_valuesE(1:end-1))+G.dx_min);
dx_valuesE = [dx_valuesE; G.dx_min];


%% Membrane
dx_valuesM = G.dx_min*ones(min(G.num_boundary_layer_dx, round(G.Lix/(2*G.dx_min))), 1);

x = min(G.num_boundary_layer_dx, round(G.Lm/(2*G.dx_min)))*G.dx_min + G.dx_min/2;
counter = 0;
current_dx = G.dx_min*increase_factor;
while x < G.Lm/2
    dx_valuesM = [dx_valuesM; current_dx];
    x = x + current_dx;
    counter = counter + 1;
    if counter == G.num_per_dx_increase
        if current_dx*increase_factor <= G.dx_max
            current_dx = current_dx*increase_factor;
        end
        counter = 0;
    end
end

if 2*sum(dx_valuesM) - G.Lm == dx_valuesM(end)
    dx_valuesM = [dx_valuesM(1:end-1); flip(dx_valuesM)];
elseif 2*sum(dx_valuesM) - G.Lm > 0
    if G.Lm-sum([dx_valuesM(1:end-1); dx_valuesM(1:end-1)]) > 0
        dx_valuesM = [dx_valuesM(1:end-1); G.Lm-sum([dx_valuesM(1:end-1); ...
            dx_valuesM(1:end-1)]); flip(dx_valuesM(1:end-1))];
    elseif G.Lm-sum([dx_valuesM(1:end-2); dx_valuesM(1:end-2)]) > 0
        dx_valuesM = [dx_valuesM(1:end-2); G.Lm-sum([dx_valuesM(1:end-2); ...
            dx_valuesM(1:end-2)]); flip(dx_valuesM(1:end-2))];
    else
        error('Not able to generate adaptive mesh. Please update the code.')
    end
elseif 2*sum(dx_valuesM) - G.Lm == 0
    dx_valuesM = [dx_valuesM; flip(dx_valuesM)];
elseif 2*sum(dx_valuesM) - G.Lm < 0
    dx_valuesM = [dx_valuesM; G.Lm-2*sum(dx_valuesM); flip(dx_valuesM)];
end

%% Active area
dx_valuesA1 = G.dx_AMPA*ones(floor((G.L_AMPA/2-G.Lsv/2-G.w_ch-4-G.Lm)/G.dx_AMPA), 1);
dx_valuesA1(end) = (G.L_AMPA/2-G.Lsv/2-G.w_ch-4-G.Lm) - sum(dx_valuesA1(1:end-1));
dx_valuesAc = G.dx_ch*ones(floor((G.w_ch)/G.dx_ch), 1);
dx_valuesAm = G.dx_min*ones(floor((G.Lm)/G.dx_min), 1);
dx_valuesA2 = (G.Lsv/2-G.w_ch/2);
dx_valuesA = [dx_valuesA1; dx_valuesAc; dx_valuesAc; dx_valuesAm; ...
    dx_valuesA2; dx_valuesAc; flip(dx_valuesA2); dx_valuesAm; ...
    dx_valuesAc; dx_valuesAc; flip(dx_valuesA1)];


%% Intracellular space
i_length = (G.Lix-G.L_AMPA)/2 - G.dx_min;
dx_valuesI = [G.dx_min; (i_length/2); (i_length/2)];


%% Collect all
dx_values = [flip(dx_valuesE); dx_valuesM; dx_valuesI; dx_valuesA; ...
    flip(dx_valuesI); dx_valuesM; dx_valuesE];
Nx = length(dx_values);

end
