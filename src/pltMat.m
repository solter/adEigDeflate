##-*- texinfo -*-
##@deftypefn{Function File}{[] = } pltMat(@var{A})
##
##Given a matrix @var{A}, plots the size of its
##entries on a logarithmic scale
##@end deftypefn

function [] = pltMat(H)
    imagesc(logNAN10(abs(H)));
    daspect([1 1 1]);
    h = colorbar;
    ytick = get(h, "ytick");
    set (h, "yticklabel", sprintf ('10^{%g}|', ytick));
end#function

function [out] = logNAN10(in)
  out = log10(in);
  out(isinf(out)) = nan; 
  out(out < -16) = nan;
endfunction
