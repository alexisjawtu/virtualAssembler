function [y] = WE_basis (x, n_vert)
  %% returns the evaluation of the four or five functions in the basis of W(E)
  %% at the point x. x should be a column vector.
  %% n_vert is to know which type of element we have
  %%
  %% TODO, tomo el promedio, despues veo si no va
  if n_vert     == 4
    p = [eye(3), x];
  elseif n_vert == 5
    p = [eye(3), x];
  elseif n_vert == 6  
    p = [eye(3), [x(1:2);0], [0;0;x(3)]];
  else
    error ('invalid number of vertices');
  end  
  y = p;
endfunction
