%% assmbl_prism: function description
function [outputs] = assmbl_prism(face_indices,face_types,Mdistances,det_M_Element,vol_pts,face_pts,P0,normalFacesE,measFacesE)
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

  clear('Mdistances');

  face_quad_coef = repmat([1;1;1;1;4;4;4;4;16],1,5);
  quadrilat = find(face_types==4);
  face_quad_coef(:,quadrilat) = face_quad_coef(:,quadrilat)*diag(measFacesE(quadrilat)/36);
  triangles = find(face_types == 3);
  face_quad_coef(:,face_indices(triangles)) = repmat([measFacesE(triangles)/3,zeros(1,6)].',1,3);

prismatic_face_points =  SEGUIR ACA






  measE                       = abs(det_M_Element)/3;

  quad_nrmlztn                = measE/2;    

  outputs = 1;
endfunction

%prismatic_face_points          = zeros(3,9,5);
prismatic_face_points(:,:,1)   = 
[0 .5 1 0 .5  1 0 .5 1; 
 0 0  0  0  0  0 0 0 0; 
 0 0  0 .5 .5 .5 1 1 1];

prismatic_face_points(:,:,2)   = 
[1 .5 0 1 .5 0 1 .5 0; 
 0 .5 1 0 .5 1 0 .5 1; 
 0 0 0 .5 .5 .5 1 1 1];

prismatic_face_points(:,:,3)   = 
[0  0 0  0  0  0 0  0 0; 
 0 .5 1  0 .5  1 0 .5 1; 
 0  0 0 .5 .5 .5 1  1 1];

prismatic_face_points(:,1:3,4) = 
[.5 .5  0; 
  0 .5 .5; 
  0  0  0];

prismatic_face_points(:,1:3,5) = 
[.5 .5  0; 
  0 .5 .5; 
  1  1  1];
