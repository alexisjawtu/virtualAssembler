%    deter         = abs(det(M_Element));
%    % TODO: with formula base x height / 3 ??;
%    measE         = meas_Refs{n_VERT} * deter;
%
%    quad_nrmlztn  = measE * quad_aux_const{n_VERT};
%    
%    centroid = mean(vertices(:,elements(el,2:n_VERT+1)),2);
%
%    diameter = max(norm(Mdistances,2,'columns'));
%
%    rescale_factor(5) = 1/diameter;
%    
%    %% obs: matrix product needs to be reshaped after:
%    face_pts    = reshape(M_Element * points_of_faces{cell_index{n_VERT}} + P(:,1), 3, max(n_face_pts{n_VERT}), nFacesE);
%    w           = zeros(3, dim_WE{n_VERT}, size(vol_pts,2)); %% polynomials of W(E)
%    potentials  = zeros(size(vol_pts,2), dim_WE{n_VERT});
%    for l = 1:size(vol_pts,2)
%      potentials(l,:) = WE_potentials(diameter,centroid,vol_pts(:,l),n_VERT);
%      w(:,:,l)        = WE_basis(diameter,centroid,vol_pts(:,l),n_VERT);  % the four w_i s.t. span {w_i:1<=i<=4} == W_E   
%    end                                                                    % zeros(3,dim_WE{n_VERT});
%    
%    int_E_w_w = zeros(dim_WE{n_VERT});   %% int_E < w_k; w_r > dx  a.k.a. G
%    for r = 1:dim_WE{n_VERT}
%      for k = 1:dim_WE{n_VERT}
%        vals            = reshape(dot(w(:,r,:), w(:,k,:), 1), 1, n_vol_pts{n_VERT});
%        int_E_w_w(r,k)  = vals*vol_wg{n_VERT};
%      end
%    end
%    
%    int_E_w_w *= quad_nrmlztn;
%
%    H         = measE;        %% before eqref(45) page 62
%    % W       = measFacesE  %% size (1 x Ndof_E). eqref(45) page 62 and page 63.
%    Hsharp    = quad_nrmlztn*potentials.'*vol_wg{n_VERT}; %% integrals of the potentials over E; int_E_Q
%                                            %% page 63 in the middle    
%    inverse_H = H\eye(size(H,1));
%    B1        = -Hsharp * inverse_H * W;  %% page 63, second half of the page.
%
%    % Boundary term. B2_{i,j} = Int_{fi} qj dS. page 62 in the middle.
%    B2        = zeros(dim_WE{n_VERT},nFacesE); 
%    % D_{j,i} = dof_j ( w_i ). See page 64 in the middle.
%    D         = zeros(nFacesE,dim_WE{n_VERT});
%    %% loop F
%    for face = 1:nFacesE
%      for pts = 1:n_face_pts{n_VERT}(face)
%        p           = WE_potentials(diameter, centroid, face_pts(:, pts, face), n_VERT);
%        w_on_faces  = WE_basis (diameter, centroid, face_pts(:, pts, face), n_VERT);
%
%        if (size(p,2) ~= dim_WE{n_VERT})
%          error('WRONG: dim of W_E and quantity of potentials are not equal.');
%          el
%          exit;
%        end
%
%        %% Next line: we are assembling the quadrature by columns
%        B2(:,face) += face_quad_coef(pts,face) * p.';
%        %% Next line: item 13 page 68. the quadrature
%        D(face,:)  += face_quad_coef(pts,face) * normalFacesE(:,face).' * w_on_faces;
%      end
%    end
%
%    B2 = B2./repmat(measFacesE,dim_WE{n_VERT},1);
%
%    if ((el == 3) || (el == 14))
%      save(strcat('element',num2str(el)),'B1','B2','measE','Hsharp','vol_pts',...
%        'centroid','diameter','inverse_H','Mdistances','int_E_w_w');
%    end
%
%    if ~all(all(size(B1)==size(B2)))
%      el
%      size(B1)
%      size(B2)
%    end
%    B            = B1 + B2;    %% see item 9 page 68

%    Pi_zero_star = int_E_w_w\B;  %% item 12, page 68. page 19 on alexis google drive.
%
%    PiW       = D*Pi_zero_star; %% item 14, page 68.
%    %% next lines: build local matrix K 11.3 page 69
%    %% computability term.
%    K_comp_a  = Pi_zero_star.' * int_E_w_w * Pi_zero_star;
%    temp_mtrx = eye(size(PiW,1)) - PiW;
%    %% stability term.
%
%    K_stab_a(5)  = rescale_factor{5} * measE * temp_mtrx.'*temp_mtrx; 
%    A = K_comp_a + K_stab_a{n_VERT};


%  test_verts = vertices(:,elements(3,2:5));
%  save('vertices_element_3', 'test_verts');
%  face1 = vertices(:,faces(10,2:4));
%  face2 = vertices(:,faces(11,2:4));
%  face3 = vertices(:,faces(12,2:4));
%  face4 = vertices(:,faces(13,2:4));
%  save('faces_element_3', 'face1', 'face2', 'face3', 'face4');
%  test_verts = vertices(:,elements(14,2:6));
%  save('vertices_element_14', 'test_verts');
%  face1 = vertices(:,faces(50,2:4));
%  face2 = vertices(:,faces(43,2:4));
%  face3 = vertices(:,faces(51,2:4));
%  face4 = vertices(:,faces(48,2:4));
%  face5 = vertices(:,faces(41,2:5));
%  save('faces_element_14', 'face1', 'face2', 'face3', 'face4', 'face5');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dim_WE            = {};  



  %meas_Refs         = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n_vol_pts         = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vol_wg            = {};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quad_aux_const    = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rescale_factor    = {};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  K_stab_a          = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  K_stab_a(4) = zeros(4);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  K_stab_a(6) = zeros(5);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rescale_factor(4) = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rescale_factor(6) = 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dim_WE(4)   = 4;  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dim_WE(5)   = 4;  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dim_WE(6)   = 5;  

  %face_pts = {};
  %face_pts(4)  = {};  
  %face_pts(5)  = {};  
  %face_pts(6)  = {};

  %meas_Refs(4) = 1/6;
  %meas_Refs(6) = 1/2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  quad_aux_const(4) = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  quad_aux_const(5) = .5;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  quad_aux_const(6) = 1/36; 

  %% quadrature points are requested only once:

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 cell_index    = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 cell_index(4) = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 cell_index(5) = 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 cell_index(6) = 3;
  