%% Copyright (C)  
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; If not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn  {Function File} {@var{A} =} assembleA (@var{filename})
%% Short Description: assembles A
%%
%% Long Description
%% 
%% shared_faces.txt reads:
%% ------------------------------
%% el1 el2 share face v1 v2 v3 v4 
%% el1 el2 share face v1 v2 v3 v4 
%% el1 el2 share face v1 v2 v3 v4  
   
%% as long as elements appear in elements_by_vertices.txt, THAT FILE 
%% determines global enumeration of faces and, at the same time,
%% local--to--global relationship of faces in elements
%%
%% faces_global.txt (3 or 4 means triangle or rectangle):
%% ------------------------------
%% 3 v1 v2 v3
%% 4 v1 v2 v3 v4 etc...

%% faces_local_to_global.txt (0 == tetra, 1 == pyr, 2 == prism)
%% ------------------------------
%% 0 f1 f2 f3 f4
%% 1 f1 f2 f3 f4 f5
%% 2 f1 f2 f3 f4 f5     (f_i is the line number in faces_global.txt)
%%
%% local traverse of faces:
%% ------------------------------
%% prism:       1: {z==0}; 2: {z==1}; 3: {y==0}; 4: {x==0}; 5: {x+y==1}
%% pyramid:     1: {z==0}; 2, 3, and 4 counterclockwise; 5: square face
%% tetrahedron: 1: {z==0}; 2: {y==0}; 3: {x==0}; 4: {x+y+z==1};

%% 
%% @seealso{assemble_[type]}
%% @end deftypefn
%% Author: Alexis Jawtuschenko.
function [res] = assembleA()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55 dict_save = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  dict_save2 = {};

  vertices              = importdata ('vertices.txt');
  elements              = importdata ('elements_by_vertices.txt');
  elements_by_faces     = importdata ('elements_by_faces.txt');
  faces                 = importdata ('faces.txt');

  vertices              = vertices.';
  num_el                = size (elements,1);
  %% num_fc is the number of faces in the whole mesh
  num_fc                = size (faces,1);

  fprintf('Global normals and measures, %d faces\n', num_fc);

  assemble_local = {};
  assemble_local(4) = @assemble_tetr;
  assemble_local(5) = @assemble_pyram;
  assemble_local(6) = @assemble_prism;
  
  global_normals       = cross(vertices(:,faces(:,3))-vertices(:,faces(:,2)),vertices(:,faces(:,4))-vertices(:,faces(:,2)));
  global_norms         = norm(global_normals,2,'columns');  
  global_normals       = global_normals./global_norms;
  global_meas_faces    = global_norms./(1 + (faces(:,1) == 3));
  to_save              = global_norms.';
  
  fprintf('Global normals and measures finished.');
  save('globalNorms.txt','to_save','num_el','num_fc');
  clear('to_save');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%dim_WE            = {};  



  %meas_Refs         = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% n_vol_pts         = {};
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vol_wg            = {};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quad_aux_const    = {};    %% adjusting the volume quadrature
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

  n_Faces     = {};
  n_Faces(4)  = 4;  
  n_Faces(5)  = 5;  
  n_Faces(6)  = 5;

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
  
  K = sparse(num_fc + num_el,num_fc + num_el);

  fprintf('Elements loop. %d\n',num_el);
  tic

  for el = 1:num_el
    h            = waitbar(el/num_el);
    n_VERT       = elements(el,1);      %% to know which type of element
    P            = vertices(:,elements(el,2:(n_VERT+1)));
    normalFacesE = global_normals(:,elements_by_faces(el,2:(n_Faces{n_VERT}+1)));
    measFacesE   = global_meas_faces(elements_by_faces(el,2:(n_Faces{n_VERT}+1)));
    faces_of_E   = elements_by_faces(el,2:(n_Faces{n_VERT}+1));
    
    face_pts = {};
    for f = 2:(n_Faces{n_VERT}+1)                   % 2:(1 + type--of--face)   
      face = vertices(:,faces(elements_by_faces(el,f),2:1+faces(elements_by_faces(el,f),1)));
      if faces(elements_by_faces(el,f),1) == 3
        face_pts(elements_by_faces(el,f)) = (face + shift(face,1,2))/2; 
      else %% see the structure of face_quad_coef in the comments in assembl_pyram
        face_pts(elements_by_faces(el,f)) = [face, (face + shift(face,1,2))/2, mean(face,2)];
      end
    end

    local = assemble_local{n_VERT}(P,faces_of_E,faces(faces_of_E,1),face_pts,normalFacesE,measFacesE);
    W     = ones(1,n_Faces{n_VERT});  %% size (1 x Ndof_E). eqref(45) page 62 and page 63.
    K(elements_by_faces(el,2:(n_Faces{n_VERT}+1)),elements_by_faces(el,2:(n_Faces{n_VERT}+1))) += local; 
    K(num_fc + el,elements_by_faces(el,2:(n_Faces{n_VERT}+1))) = W;
    K(elements_by_faces(el,2:(n_Faces{n_VERT}+1)),num_fc + el) = W.';

    clear('Mdistances');

  end    %% end of main for
  close(h);
  res = K;
  a=toc;
  toc
  csvwrite('TIME',a);

endfunction

%% T: just for plotting purposes
function [y] = T(M,xE,x)
  y = M*x+xE;
endfunction








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
%    face_pts    = reshape(M_Element * points_of_faces{cell_index{n_VERT}} + P(:,1), 3, max(n_face_pts{n_VERT}), n_Faces{n_VERT});
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
%    B2        = zeros(dim_WE{n_VERT},n_Faces{n_VERT}); 
%    % D_{j,i} = dof_j ( w_i ). See page 64 in the middle.
%    D         = zeros(n_Faces{n_VERT},dim_WE{n_VERT});
%    %% loop F
%    for face = 1:n_Faces{n_VERT}
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