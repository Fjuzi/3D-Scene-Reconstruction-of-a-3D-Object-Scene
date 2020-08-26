function A = internal_parameters_solve_vmmc(H)
% internal_parameters_solve_vmmc(H) returns the intrinsic matrix following the
% Zhang's camera calibration method. H is a cell-structure where each
% field H{i} is an homography between the pattern and each image. 
%
% This method computes A through the conic B = A'*inv(A). Constraints over
% the conic and H{i} are derived into an homogeneous system of equations 
% defined by V·b=0, which can be solved through svd decomposition.
% More details at http://research.microsoft.com/en-us/um/people/zhang/Papers/TR98-71.pdf
%
%
% Returns the 3x3 matrix (the intrinsic matrix A)
%
% Juan C. SanMiguel, Universidad Autonoma de Madrid, November 2012

%computation of V for all the images
V=[];
for i=1:numel(H)
    V = [V; vfun(H{i},1,2)'; (vfun(H{i},1,1)-vfun(H{i},2,2))'];
end

%solution through svd
[ub,sb,vb] = svd(V);
b = vb(:,6);

%extraction of the intrinsic parameters according to B matrix
bn = b(2)*b(4)-b(1)*b(5);
bd = b(1)*b(3)-b(2)^2;

v0 = bn/bd;                             % v0 coordinate of the principal point.
lambda = b(6)-(b(4)^2+v0*bn)/b(1);      
alfa = sqrt(abs(lambda/b(1)));          % scale factor for the u axis.
beta = sqrt(lambda*b(1)/bd);            % scale factor for the v axis.
gamma = -b(2)*alfa^2*beta/lambda;       % parameter representing a non-perpendicular image plane.
u0 = gamma*v0/alfa-b(4)*alfa^2/lambda;  % u0 coordinate of the principal point.

A = [alfa gamma u0
     0	  beta  v0
     0	  0	    1];

function v=vfun(h,i,j)
v = [h(1,i)*h(1,j)
    h(1,i)*h(2,j)+h(2,i)*h(1,j)
    h(2,i)*h(2,j)
    h(3,i)*h(1,j)+h(1,i)*h(3,j)
    h(3,i)*h(2,j)+h(2,i)*h(3,j)
    h(3,i)*h(3,j)];
return