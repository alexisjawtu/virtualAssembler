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
%% @deftypefn  {Function File} {@var{A} =} assemble (@var{filename})
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
function [K,F,num_fc,num_el] = assemble()

  vertices              = importdata ('vertices.txt');
  elements              = importdata ('elements_by_vertices.txt');
  elements_by_faces     = importdata ('elements_by_faces.txt');
  faces                 = importdata ('faces.txt');

  vertices              = vertices.';
  num_el                = size (elements,1)
  num_fc                = size (faces,1)

  fprintf('Global normals and measures, %d faces\n', num_fc);

  assemble_local    = {};
  assemble_local(4) = @assemble_tetr;
  assemble_local(5) = @assemble_pyram;
  assemble_local(6) = @assemble_prism;
  
  global_normals    = cross(vertices(:,faces(:,3))-vertices(:,faces(:,2)),vertices(:,faces(:,4))-vertices(:,faces(:,2)));
  global_norms      = norm(global_normals,2,'columns');  
  global_normals    = global_normals./global_norms;
  global_meas_faces = global_norms./(1 + (faces(:,1) == 3));
  to_save           = global_norms.';
  
  fprintf('Global normals and measures finished.');
  save('globalNorms.txt','to_save','num_el','num_fc');
  clear('to_save');

  n_Faces     = {};
  n_Faces(4)  = 4;  
  n_Faces(5)  = 5;  
  n_Faces(6)  = 5;

  K = sparse(num_fc + num_el,num_fc + num_el);
  F = sparse(num_fc + num_el,1);

  fprintf('Elements loop. %d\n',num_el);
  tic

  for el = 1:num_el
    h            = waitbar(el/num_el);
    n_VERT       = elements(el,1);      %% to know which type of element
    nFacesE      = n_Faces{n_VERT};
    P            = vertices(:,elements(el,2:(n_VERT+1)));
    faces_of_E   = elements_by_faces(el,2:(nFacesE+1));
    normalFacesE = global_normals(:,faces_of_E);
    measFacesE   = global_meas_faces(faces_of_E);
    face_types   = faces(faces_of_E,1);
    
    face_pts = {};
    for f = 2:(nFacesE+1)                   % 2:(1 + type--of--face)   
      face = vertices(:,faces(elements_by_faces(el,f),2:1+faces(elements_by_faces(el,f),1)));
      if faces(elements_by_faces(el,f),1) == 3
        face_pts(elements_by_faces(el,f)) = (face + shift(face,1,2))/2; 
      else %% see the structure of face_quad_coef in the comments in assembl_pyram
        %% after testing
        face_pts(elements_by_faces(el,f)) = [face, (face + shift(face,1,2))/2, mean(face,2)];
      end
    end
    [local_matrix, local_F]   = assemble_local{n_VERT}(P,faces_of_E,face_types,face_pts,normalFacesE,measFacesE);
    W                         = ones(1,nFacesE);  %% size (1 x Ndof_E). eqref(45) page 62 and page 63.
    K(faces_of_E,faces_of_E) += local_matrix; 
    K(num_fc + el,faces_of_E) = W;
    K(faces_of_E,num_fc + el) = W.';
    F(num_fc + el,1)          = local_F;

  end    %% end of main for
  close(h);
  a=toc;
  toc
  csvwrite('TIME',a);
  K_full = full(K);
  F_full = full(F);
  save ('assembled.txt', 'K_full','F_full')

endfunction

%% T: just for plotting purposes
function [y] = T(M,xE,x)
  y = M*x+xE;
endfunction