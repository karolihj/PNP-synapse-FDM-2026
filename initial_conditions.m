function c0 = initial_conditions(x, y, z, G, c_idx)
%c0 = initial_conditions(x, y, z, G, c_idx)

c0 = G.c0e(c_idx)*ones(size(x));
c0(find(in_membrane(x,y,z,G))) = 0;
c0(find(in_vesicle(x,y,z,G))) = G.c0v(c_idx);
c0(find(in_cell(x,y,z,G))) = G.c0i(c_idx);

end

