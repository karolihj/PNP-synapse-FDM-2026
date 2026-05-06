function [x, y, z, x_idx, y_idx, z_idx] = x_y_z_from_idx(idx, G)
%[x, y, z, x_idx, y_idx, z_idx] = x_y_z_from_idx(idx, G)

x_idx = rem((idx-1), G.Nx)+1;
y_idx = rem(floor((idx-1)/G.Nx), G.Ny)+1;
z_idx = floor((idx-1)/(G.Nx*G.Ny))+1;

x = zeros(length(idx), 1);
y = zeros(length(idx), 1);
z = zeros(length(idx), 1);
for n=1:length(idx)
    x(n) = sum(G.dx_values(1:x_idx(n)-1)) + G.dx_values(x_idx(n))/2;
    y(n) = sum(G.dy_values(1:y_idx(n)-1)) + G.dy_values(y_idx(n))/2;
    z(n) = sum(G.dz_values(1:z_idx(n)-1)) + G.dz_values(z_idx(n))/2;
end

end

