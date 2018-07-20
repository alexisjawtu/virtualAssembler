%% assmble_prism: function description
function [outputs] = assmble_prism(vertices,face_indices,face_types,face_pts,normalFacesE,measFacesE)
  %% for now, Vh(prism) == Wh(prism)
  %
  % TODO: make a benchmark for the initialization of pts_of_faces each time
  % vs. passing it to the function
  %
  % Explanation of vol_pts: if E were the "reference prism", we would have
  %
  % vol_pts = [.5  0 0; .5 .5 0; 0 .5 0; .5 0 .5; .5 .5 .5; 0 .5 .5; .5 0 1; .5 .5 1; 0 .5 1];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %%  this was in assembleA.m. It will be useful in future versions of the program.

  %
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  n_vertices   = 6;
  n_vol_pts    = 9;
  dim_Vh       = 5;
  n_vertices   = 6;

  vol_weights  = [1; 4; 1; 1; 4; 1; 1; 4; 1];

  we_basis      = zeros(3, dim_Vh, n_vol_pts);
  M_Element     = [vertices(:,2)-vertices(:,1),vertices(:,3)-vertices(:,1),vertices(:,4)-vertices(:,1)];
  det_M_Element = det(M_Element);
  measE         = abs(det_M_Element)/2;
  quad_nrmlztn  = measE/36;    

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
  
  int_E_w_w = zeros(dim_Vh);   %% int_E < w_k; w_r > dx
  for r = 1:dim_Vh
    for k = 1:dim_Vh
      %% ( (wr*wk)(p1), ..., (wr*wk)(p9) )
      vals           = reshape(dot(we_basis(:,r,:),we_basis(:,k,:),1),1,n_vol_pts);
      int_E_w_w(r,k) = vals*vol_weights;
    end
  end
  
  int_E_w_w *= quad_nrmlztn; 
  outputs    = int_E_w_w;
endfunction