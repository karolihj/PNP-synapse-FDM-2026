function [in_area] = in_vesicle(x,y,z, G)
%[in_area] = in_vesicle(x,y,z,G)

in_area = (x > G.Lex+G.Lm+(G.Lix-G.Lsv)/2).*(x < G.Lex+G.Lm+(G.Lix+G.Lsv)/2)...
    .*(z > G.Lez+G.Lm+(G.Liz-G.Lsv)/2).*(z < G.Lez+G.Lm+(G.Liz+G.Lsv)/2)...
    .*(y > G.Liy+G.L_cleft+2*G.Lm).*(y < G.Liy+G.L_cleft+2*G.Lm+G.Lsv);

end

