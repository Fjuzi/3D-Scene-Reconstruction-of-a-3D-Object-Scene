function [R T] = external_parameters_solve_vmmc(A, H)
% [R T]=external_parameters_solve_vmmc(A, H) returns the extrinsic matrix [R, T] 
% following the Zhang's camera calibration method.  More details at 
% http://research.microsoft.com/en-us/um/people/zhang/Papers/TR98-71.pdf
% 
% Input:
%   - A is a 3x3 matrix defining the intrinsic camera parameters
%   - H is a cell-structure where each field H{i} is 3x3 matrix which describes 
%     an homography between the pattern and each image. 
%
% Output:
%   - R is a 3x3 matrix representing the rotation parameters
%   - T is a 3x1
%
% Juan C. SanMiguel, Universidad Autonoma de Madrid, November 2012

inva = pinv(A);                     % inverse of intrinsic parameter matrix
x0 = [];
% 
for i=1:numel(H)
    h = H{i};                      % homography H (matriz de proyección).
    
    %Rotation matrix estimation
    lamb = 1/norm(inva*h(:,1));
    r1 = lamb*inva*h(:,1);          % rotaton vector r1.
    r2 = lamb*inva*h(:,2);          % rotaton vector r2.
    r3 = cross(r1,r2);              % rotaton vector r3.
    
    %Orthogonality Enforcement    
    R_temp = [r1 r2 r3];
    [u,s,v] = svd(R_temp);
    Q{i} = u*v';
    
    %computation of the rotation matrix
    ang{i} = matrot_vmmc(Q{i});    
    R{i} = matrot_vmmc(ang{i});            
    R{i} = R{i}(1:3,1:3);
    
    % Translation Vector Estimation
    T{i} = lamb*inva*h(:,3);        % translation vector t.    
end