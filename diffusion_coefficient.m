function D = diffusion_coefficient(x, y, z, G, c_idx)
%D = diffusion_coefficient(x, y, z, G, c_idx)

D = G.D(c_idx)*ones(size(x));
D(find(in_membrane(x,y,z,G))) = 0;
D(find(in_cleft(x,y,z,G))) = G.D(c_idx)/G.kappa_cleft(c_idx);

if isfield(G, 'open_vesicle') && G.open_vesicle
    D(find(in_vesicle_opening(x,y,z,G))) = G.D(c_idx)/G.kappa_cleft(c_idx);
end

end

