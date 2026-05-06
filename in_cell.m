function [in_area] = in_cell(x,y,z,G)
%[in_area] = in_cell(x,y,z,G)

in_lower = (x > G.Lex + G.Lm).*(x < G.Lex + G.Lm + G.Lix).*...
    (z > G.Lez + G.Lm).*(z < G.Lez + G.Lm + G.Liz).*(y < G.Liy);

in_upper = (x > G.Lex + G.Lm).*(x < G.Lex + G.Lm + G.Lix).*...
    (z > G.Lez + G.Lm).*(z < G.Lez + G.Lm + G.Liz).*...
    (y > G.Liy+2*G.Lm+G.L_cleft);

in_v = in_vesicle(x,y,z,G);
in_m = in_membrane(x,y,z,G);

in_upper(in_v == 1) = 0;
in_upper(in_m == 1) = 0;

in_area = in_lower + in_upper;

end