function [y] = WE_potentials(x, n_vert)  
  if n_vert == 4                 %  obs: dot is faster
    p = [(x).', dot(x,x)/(2)];
%   p = [(x).', diameter*0.5*sum(((x)/diameter).^2)];
  elseif n_vert == 5
    p = [(x).', dot(x,x)/(2)];
  elseif n_vert == 6
    p = [(x).', dot((x)([1:2]),(x)([1:2]))/(2), ((x)(3))^2/(2)];
  else
    error ('invalid number of vertices');
  end
  y = p;
endfunction
