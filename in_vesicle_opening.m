function [in_area] = in_vesicle_opening(x,y,z,G)
%[in_area] = in_vesicle_opening(x,y,z)

in_area = (x > G.Lex+G.Lm+(G.Lix-G.w_vo)/2).*(x < G.Lex+G.Lm+(G.Lix+G.w_vo)/2)...
    .*(z > G.Lez+G.Lm+(G.Liz-G.w_vo)/2).*(z < G.Lez+G.Lm+(G.Liz+G.w_vo)/2)...
    .*(y >= G.Liy+G.L_cleft+G.Lm).*(y <= G.Liy+G.L_cleft+2*G.Lm);

end