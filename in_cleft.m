function [in_area] = in_cleft(x,y,z, G)
%[in_area] = in_cleft(x,y,z,G)

in_area = (x >= G.Lex).*(x <= G.Lex+G.Lix+2*G.Lm).*(z >= G.Lez).*...
    (z <= G.Lez+G.Liz+2*G.Lm).*(y > G.Liy+G.Lm).*(y < G.Liy+G.L_cleft+G.Lm);

end