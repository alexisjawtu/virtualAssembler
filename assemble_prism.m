%% assmbl_prism: function description
function [outputs] = assmbl_prism(vertices,face_indices,face_types,face_pts,normalFacesE,measFacesE)
  %% for now, Vh(prism) == Wh(prism)
  %
  %% measFacesE in R(1x5)
  %
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
  
  %% this was in assembleA.m. It will be useful in future versions of the program.
  %      face_pts = {};
  %
  %      for f = 2:(n_Faces{n_VERT}+1)                   
  %      % list of vertices of the face (3x3 or 3x4)       2:(1 + type--of--face)   
  %        face = vertices(:,faces(elements_by_faces(el,f),2:(1+faces(elements_by_faces(el,f),1))));
  %        if faces(elements_by_faces(el,f),1) == 3
  %          % TODO: use the mean instead of (a+b)/2
  %          face_pts(elements_by_faces(el,f)) = (face + shift(face,1,2))/2; 
  %        else %% see the structure of face_quad_coef in the comments in assembl_pyram
  %          face_pts(elements_by_faces(el,f)) = [face, (face + shift(face,1,2))/2, mean(face,2)];
  %        end
  %      end

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
			**ver todo assembleA.m y borrar todo lo que no se use y controlar cada linea
			** ver tema face_pts en programa prismas

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