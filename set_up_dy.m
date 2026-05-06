function [dy_values, Ny] = set_up_dy(G, channels)
% [dy_values, Ny] = set_up_dy(G)
% Set up adaptive mesh in the y-direction

%% From bottom to lower potassium channels
dy_values1 = [((channels{3}.channel_start_y-G.dy_min)/2)*ones(2,1); G.dy_min];


%% From potassium channels to sodium channels
l1 = channels{7}.channel_start_y - (channels{3}.channel_start_y + channels{3}.ly);
dy_values2 = (l1/2)*ones(2,1);


%% From sodium channels to membrane
l1 = G.Liy - (channels{7}.channel_start_y + channels{7}.ly);
y =  G.num_boundary_layer_dy*G.dy_min + G.dy_min/2;
dy_values3 = G.dy_min*ones(G.num_boundary_layer_dy, 1);
current_dy = G.dy_min*G.increase_factor_y;
counter = 0;
while y < l1 
    dy_values3 = [dy_values3; current_dy];
    y = y + current_dy;
    counter = counter + 1;
    if counter == G.num_per_dy_increase
        if current_dy*G.increase_factor_y <= G.dy_max
            current_dy = current_dy*G.increase_factor_y;
        end
        counter = 0;
    end
end

dy_values3(end) = l1 -sum(dy_values3(1:end-1));


%% From sodium channels to membrane presynaptic cell
% Synaptic vesicle
l1 = G.Lsv - 2*G.num_boundary_layer_dy*G.dy_min;
dy_values3a = G.dy_min*G.increase_factor_y;
current_dy = dy_values3a*G.increase_factor_y;
y = dy_values3a/2;

while y < l1/2
    dy_values3a = [dy_values3a; current_dy];
    y = y + current_dy;
    counter = counter + 1;
    if counter == G.num_per_dy_increase
        if current_dy*G.increase_factor_y <= G.dy_max
            current_dy = current_dy*G.increase_factor_y;
        end
        counter = 0;
    end
end

dy_values3a = dy_values3a(1:end-1);
dy_values3a(end) = l1/2 -sum(dy_values3a(1:end-1));
dy_values3a = [G.dy_min*ones(G.num_boundary_layer_dy, 1); dy_values3a; ...
    flip(dy_values3a); G.dy_min*G.num_boundary_layer_dy];


% After synaptic vesicle
l1 = G.Liy - (channels{7}.channel_start_y + channels{7}.ly) - G.Lsv + G.Lm;
dy_values3b = G.dy_min*G.num_boundary_layer_dy;
current_dy = dy_values3b*G.increase_factor_y*2;
y =  dy_values3b/2;
counter = 0;
while y < l1 
    dy_values3b = [dy_values3b; current_dy];
    y = y + current_dy;
    counter = counter + 1;
    if counter == G.num_per_dy_increase
        if current_dy*G.increase_factor_y <= G.dy_max
            current_dy = current_dy*G.increase_factor_y;
        end
        counter = 0;
    end
end

dy_values3b(end) = l1 -sum(dy_values3b(1:end-1));

%% Cleft
dy_valuesI = G.dy_min*ones(min(G.num_boundary_layer_dy, round(G.L_cleft/(2*G.dy_min))), 1);

y = min(G.num_boundary_layer_dy, round(G.L_cleft/(2*G.dy_min)))*G.dy_min + G.dy_min/2;
counter = 0;
current_dy = G.dy_min*G.increase_factor_y;
while y < G.L_cleft/2
    dy_valuesI = [dy_valuesI; current_dy];
    y = y + current_dy;
    counter = counter + 1;
    if counter == G.num_per_dy_increase
        if current_dy*G.increase_factor_y <= G.dy_max
            current_dy = current_dy*G.increase_factor_y;
        end
        counter = 0;
    end
end

if 2*sum(dy_valuesI) - G.L_cleft == dy_valuesI(end)
    dy_valuesI = [dy_valuesI(1:end-1); flip(dy_valuesI)];
elseif 2*sum(dy_valuesI) - G.L_cleft > 0
    if G.L_cleft-sum([dy_valuesI(1:end-1); dy_valuesI(1:end-1)]) > 0
        dy_valuesI = [dy_valuesI(1:end-1); G.L_cleft-sum([dy_valuesI(1:end-1); ...
            dy_valuesI(1:end-1)]); flip(dy_valuesI(1:end-1))];
    else
        error('Not able to generate adaptive mesh. Please update the code.')
    end
elseif 2*sum(dy_valuesI) - G.L_cleft == 0
    dy_valuesI = [dy_valuesI; flip(dy_valuesI)];
elseif 2*sum(dy_valuesI) - G.L_cleft < 0
    dy_valuesI = [dy_valuesI; G.L_cleft-2*sum(dy_valuesI); flip(dy_valuesI)];
end


%% Membrane
dy_valuesM = G.dy_min*ones(min(G.num_boundary_layer_dy, round(G.Liy/(2*G.dy_min))), 1);

y = min(G.num_boundary_layer_dy, round(G.Lm/(2*G.dy_min)))*G.dy_min + G.dy_min/2;
counter = 0;
current_dy = G.dy_min*G.increase_factor_y;
while y < G.Lm/2
    dy_valuesM = [dy_valuesM; current_dy];
    y = y + current_dy;
    counter = counter + 1;
    if counter == G.num_per_dy_increase
        if current_dy*G.increase_factor_y <= G.dy_max
            current_dy = current_dy*G.increase_factor_y;
        end
        counter = 0;
    end
end

if 2*sum(dy_valuesM) - G.Lm == dy_valuesM(end)
    dy_valuesM = [dy_valuesM(1:end-1); flip(dy_valuesM)];
elseif 2*sum(dy_valuesM) - G.Lm > 0
    if G.Lm-sum([dy_valuesM(1:end-1); dy_valuesM(1:end-1)]) > 0
        dy_valuesM = [dy_valuesM(1:end-1); G.Lm-sum([dy_valuesM(1:end-1); ...
            dy_valuesM(1:end-1)]); flip(dy_valuesM(1:end-1))];
    elseif G.Lm-sum([dy_valuesM(1:end-2); dy_valuesM(1:end-2)]) > 0
        dy_valuesM = [dy_valuesM(1:end-2); G.Lm-sum([dy_valuesM(1:end-2); ...
            dy_valuesM(1:end-2)]); flip(dy_valuesM(1:end-2))];
    else
        error('Not able to generate adaptive mesh. Please update the code.')
    end
elseif 2*sum(dy_valuesM) - G.Lm == 0
    dy_valuesM = [dy_valuesM; flip(dy_valuesM)];
elseif 2*sum(dy_valuesM) - G.Lm < 0
    dy_valuesM = [dy_valuesM; G.Lm-2*sum(dy_valuesM); flip(dy_valuesM)];
end


%% Collect all
dy_values = [flip(dy_values1); G.dy_ch*ones(round(channels{3}.ly/G.dy_ch),1); ...
    dy_values2; G.dy_ch*ones(round(channels{7}.ly/G.dy_ch),1); flip(dy_values3); ...
    dy_valuesM; dy_valuesI; dy_valuesM; 
    dy_values3a; G.Lm; dy_values3b;
    G.dy_ch*ones(round(channels{7}.ly/G.dy_ch),1); flip(dy_values2); ...
    G.dy_ch*ones(round(channels{3}.ly/G.dy_ch),1); dy_values1];
Ny = length(dy_values);


end