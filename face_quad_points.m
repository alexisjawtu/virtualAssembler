%% face_quad_points: function description
function [p] = face_quad_points()

  %% faces are ordered according the txts and the reference elements.

  tetrahedral_face_points        = zeros(3,3,4);
  tetrahedral_face_points(:,:,1) = [.5 .5 0; 0 .5 .5; 0 0 0];
  tetrahedral_face_points(:,:,2) = [0 0 0; .5 .5 0; 0 .5 .5];
  tetrahedral_face_points(:,:,3) = [.5 .5 0; 0 0 0; 0 .5 .5];
  tetrahedral_face_points(:,:,4) = [.5 .5 0; 0 .5 .5; .5 0 .5];

  pyramidal_face_points          = zeros(3,9,5);
  pyramidal_face_points(:,:,1)   = [0 .5 1 0 .5 1 0 .5 1; 0 0 0 .5 .5 .5 1 1 1; 0 0 0 0 0 0 0 0 0];
  pyramidal_face_points(:,1:3,2) = [0 0 0; .5 .5 0; 0 .5 .5];

  pyramidal_face_points(:,1:3,3) = [.5 .5 0; 1 .5 .5; 0 .5 .5];
  pyramidal_face_points(:,1:3,4) = [1 .5 .5; .5 .5 0; 0 .5 .5];
  pyramidal_face_points(:,1:3,5) = [.5 .5 0; 0 0 0; 0 .5 .5];

  prismatic_face_points          = zeros(3,9,5);
  prismatic_face_points(:,:,1)   = [0 .5 1 0 .5  1 0 .5 1; 0 0  0  0  0  0 0 0 0; 0 0  0 .5 .5 .5 1 1 1];
  prismatic_face_points(:,:,2)   = [1 .5 0 1 .5 0 1 .5 0; 0 .5 1 0 .5 1 0 .5 1; 0 0 0 .5 .5 .5 1 1 1];
  prismatic_face_points(:,:,3)   = [0 0 0 0 0 0 0 0 0; 0 .5 1 0 .5 1 0 .5 1; 0 0 0 .5 .5 .5 1 1 1];
  prismatic_face_points(:,1:3,4) = [.5 .5 0; 0 .5 .5; 0 0 0];
  prismatic_face_points(:,1:3,5) = [.5 .5 0; 0 .5 .5; 1 1 1];

  p = {
    tetrahedral_face_points,
    pyramidal_face_points,
    prismatic_face_points
  };
endfunction