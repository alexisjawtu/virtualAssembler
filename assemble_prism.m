%% assmbl_prism: function description
function [outputs] = assmbl_prism(vertices,face_indices,face_types,face_pts,normalFacesE,measFacesE)
  %% for now, Vh(prism) == Wh(prism)
  %
  %% measFacesE in R(1x5)
  %
  %% face_pts is a cell. face_pts comes with five GLOBAL UNIQUE face indices,
  %% those of elements_by_faces.txt
  %
  % TODO: hacer un benchmark de inicializar pts_of_faces cada vez vs. pasarlo a la funcion
  %
  %%   face_quad_coef =    |meas(fi)/3| in case of triangle
  %                        |meas(fi)/3|
  %                        |meas(fi)/3|
  %
  %%   face_quad_coef =    |1 |*(1/36) in case of rectangle
  %                        |1 |
  %                        |1 |
  %                        |1 |
  %                        |4 |
  %                        |4 |
  %                        |4 |
  %                        |4 |
  %                        |16|	

  %  face_quad_coef              = repmat([1;1;1;1;4;4;4;4;16],1,5);
  %  quadrilat                   = find(face_types == 4);
  %  triangles                   = find(face_types == 3);    
  %% TODO: instead of measFacesE(triangles) do: measFacesE(map(triangles)) 
  %% so that we don't care of "preserving orders"
  %  face_quad_coef(:,quadrilat) = face_quad_coef(:,quadrilat)*diag(measFacesE(quadrilat)/36);
  %  face_quad_coef(:,triangles) = [repmat(measFacesE(triangles)/3,3,1);zeros(6,2)];

  n_vertices   = 6;
  n_vol_pts    = 9;
  dim_Vh       = 5;
  n_vertices   = 6;
  vol_weights  = ones(n_vol_pts,1);
  we_basis     = zeros(3, dim_Vh, n_vol_pts);
  M_Element   = [vertices(:,2)-vertices(:,1),vertices(:,3)-vertices(:,1),vertices(:,4)-vertices(:,1)];
  det_M_Element = det(M_Element);
  measE        = abs(det_M_Element)/3;
  quad_nrmlztn = measE/2;    

  %% TODO: benchmark between this and (P1 + P2)/2
  vol_pts      = zeros(3,n_vol_pts);
  vol_pts(:,1) = mean([vertices(:,1),vertices(:,3)],2);
  vol_pts(:,2) = mean([vertices(:,2),vertices(:,3)],2);
  vol_pts(:,3) = mean([vertices(:,2),vertices(:,1)],2);
  vol_pts(:,7) = mean([vertices(:,4),vertices(:,6)],2);
  vol_pts(:,8) = mean([vertices(:,5),vertices(:,6)],2);
  vol_pts(:,9) = mean([vertices(:,5),vertices(:,4)],2);
  vol_pts(:,4) = mean([vol_pts(:,7),vol_pts(:,1)],2);
  vol_pts(:,5) = mean([vol_pts(:,8),vol_pts(:,2)],2);
  vol_pts(:,6) = mean([vol_pts(:,9),vol_pts(:,3)],2);

  for l = 1:n_vol_pts
    we_basis(:,:,l) = WE_basis(vol_pts(:,l),n_vertices);   
  end                                                                    
  
SEGUIR ACA: 
			**ver si la construccion de caras se puede hacer dentro de las
			funciones elementales
			**ver todo assembleA.m y borrar todo lo que no se use y controlar cada linea

  int_E_w_w = zeros(dim_Vh);   %% int_E < w_k; w_r > dx
  for r = 1:dim_Vh
    for k = 1:dim_Vh
      %% ( (wr*wk)(p1), ..., (wr*wk)(p9) )
      vals            = reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts);
      int_E_w_w(r,k)  = vals*vol_weights;
    end
  end
  
  int_E_w_w *= quad_nrmlztn; 
  outputs    = int_E_w_w;
endfunction