function [dz_values, Nz] = set_up_dz(G)
% [dz_values, Nz] = set_up_dz(G)
% Set up adaptive mesh in the z-direction

increase_factor = G.increase_factor_z;

%% Extracellular space
dz_valuesE = G.dz_min*ones(G.num_boundary_layer_dz, 1);
z = G.num_boundary_layer_dz*G.dz_min + G.dz_min/2;
counter = 0;
if G.dz_max >= G.dz_min*increase_factor
    current_dz = G.dz_min*increase_factor;
else
    current_dz = G.dz_min;
end
while z < G.Lez -G.dz_min
    dz_valuesE = [dz_valuesE; current_dz];
    z = z + current_dz;
    counter = counter + 1;
    if counter == G.num_per_dz_increase
        if current_dz*increase_factor <= G.dz_max
            current_dz = current_dz*increase_factor;
        end
        counter = 0;
    end
end

dz_valuesE(end) = G.Lez -(sum(dz_valuesE(1:end-1))+G.dz_min);
dz_valuesE = [dz_valuesE; G.dz_min];


%% Membrane
dz_valuesM = G.dz_min*ones(min(G.num_boundary_layer_dz, round(G.Liz/(2*G.dz_min))), 1);

z = min(G.num_boundary_layer_dz, round(G.Lm/(2*G.dz_min)))*G.dz_min + G.dz_min/2;
counter = 0;
current_dz = G.dz_min*increase_factor;
while z < G.Lm/2
    dz_valuesM = [dz_valuesM; current_dz];
    z = z + current_dz;
    counter = counter + 1;
    if counter == G.num_per_dz_increase
        if current_dz*increase_factor <= G.dz_max
            current_dz = current_dz*increase_factor;
        end
        counter = 0;
    end
end

if 2*sum(dz_valuesM) - G.Lm == dz_valuesM(end)
    dz_valuesM = [dz_valuesM(1:end-1); flip(dz_valuesM)];
elseif 2*sum(dz_valuesM) - G.Lm > 0
    if G.Lm-sum([dz_valuesM(1:end-1); dz_valuesM(1:end-1)]) > 0
        dz_valuesM = [dz_valuesM(1:end-1); G.Lm-sum([dz_valuesM(1:end-1); ...
            dz_valuesM(1:end-1)]); flip(dz_valuesM(1:end-1))];
    elseif G.Lm-sum([dz_valuesM(1:end-2); dz_valuesM(1:end-2)]) > 0
        dz_valuesM = [dz_valuesM(1:end-2); G.Lm-sum([dz_valuesM(1:end-2); ...
            dz_valuesM(1:end-2)]); flip(dz_valuesM(1:end-2))];
    else
        error('Not able to generate adaptive mesh. Please update the code.')
    end
elseif 2*sum(dz_valuesM) - G.Lm == 0
    dz_valuesM = [dz_valuesM; flip(dz_valuesM)];
elseif 2*sum(dz_valuesM) - G.Lm < 0
    dz_valuesM = [dz_valuesM; G.Lm-2*sum(dz_valuesM); flip(dz_valuesM)];
end

%% Active area
dz_valuesA1 = G.dz_AMPA*ones(floor((G.L_AMPA/2-G.Lsv/2-G.w_ch-4-G.Lm)/G.dz_AMPA), 1);
dz_valuesA1(end) = (G.L_AMPA/2-G.Lsv/2-G.w_ch-4-G.Lm) - sum(dz_valuesA1(1:end-1));
dz_valuesAc = G.dz_ch*ones(floor((G.w_ch)/G.dz_ch), 1);
dz_valuesAm = G.dz_min*ones(floor((G.Lm)/G.dz_min), 1);
dz_valuesA2 = (G.Lsv/2-G.w_ch/2);
dz_valuesA = [dz_valuesA1; dz_valuesAc; dz_valuesAc; dz_valuesAm; ...
    dz_valuesA2; dz_valuesAc; flip(dz_valuesA2); dz_valuesAm; ...
    dz_valuesAc; dz_valuesAc; flip(dz_valuesA1)];

%% Intracellular space
i_length = (G.Liz-G.L_AMPA)/2 -G.dz_min;
dz_valuesI = [G.dz_min; (i_length/2); (i_length/2)];

%% Collect all
dz_values = [flip(dz_valuesE); dz_valuesM; dz_valuesI; dz_valuesA; ...
    flip(dz_valuesI); dz_valuesM; dz_valuesE];
Nz = length(dz_values);

end