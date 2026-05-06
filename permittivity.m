function eps = permittivity(x, y, z, G)
%eps = permittivity(x, y, z, G)

on_membrane = in_membrane(x,y,z,G);
eps = G.eps*(~on_membrane) + G.eps_mem*on_membrane;

end

