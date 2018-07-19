%% assmbl_tetr: function description
function [outputs] = assmbl_tetr(vertices,face_indices,face_types,face_pts,normalFacesE,measFacesE)
  %% for now, Vh(tetra) == Wh(tetra)
  %
  %% measFacesE in R(1x5)
  %
  %  TODO: hacer un benchmark de inicializar pts_of_faces cada vez vs. 
  %  pasarlo a la funcion.
  %
  %% face_quad_coef =    |meas(fi)/3| in case of triangle
  %                      |meas(fi)/3|
  %                      |meas(fi)/3|
  % face_quad_coef  = repmat(measFacesE/3,3,1);
  dim_Wh          = 4; %% Wh == Vh in this case
  n_vertices      = 4;
  n_vol_pts       = 4;
  vol_pts = zeros(3,n_vol_pts); 
  M_Element       = [vertices(:,2)-vertices(:,1), vertices(:,3)-vertices(:,1), vertices(:,4)-vertices(:,1)];
  quad_nrmlztn    = abs(det(M_Element))/24; %% <--- measE/4;
  we_basis        = zeros(3, dim_Wh, n_vol_pts);
  % coefficients tetrahedral cubature on 4 points Stroud 1971 Chap 8.8
  const_a         = .58541019662496852;
  const_b         = .1381966011250105;
  
  % tetrahedral cubature points from Stroud 1971 Chap 8.8
  % direct from the physical vertices      
  vol_pts(:,1) = const_a*vertices(:,1) + const_b*sum(vertices(:,[3,4,5]),2);
  vol_pts(:,2) = const_a*vertices(:,2) + const_b*sum(vertices(:,[2,4,5]),2);
  vol_pts(:,3) = const_a*vertices(:,3) + const_b*sum(vertices(:,[2,3,5]),2);
  vol_pts(:,4) = const_a*vertices(:,4) + const_b*sum(vertices(:,[2,3,4]),2);

  for l = 1:n_vol_pts
    we_basis(:,:,l)    = WE_basis(vol_pts(:,l),n_vertices);   
  end                                                                    
  
  int_E_w_w = zeros(dim_Wh);   %% int_E < w_k; w_r > dx
  for r = 1:dim_Wh
    for k = 1:dim_Wh
      %% |E|/4 * sum ( (wr*wk)(p1), ..., (wr*wk)(p4) )
      int_E_w_w(r,k) = sum(reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts));
    end
  end
 
  outputs   = int_E_w_w * quad_nrmlztn;
endfunction
