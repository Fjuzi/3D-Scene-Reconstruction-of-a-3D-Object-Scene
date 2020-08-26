function xy = get_user_points_vmmc(ima)
% GET USER PONITS returns a list of ppoints selected by the user
%   Points should be selected by clicking the LEFT mouse button.
%   The last point should be selected with the RIGHT mouse button.
%
%   Returns a 2xN matrix of output vectors
%
% Jesus Bescos, Universidad Autonoma de Madrid, November 2012

h=figure('Name','Get User Points');
imshow(ima);
hold on;
xy = [];
np = 0;
% Loop, picking up the points.
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'r*')
    np = np+1;
    xy(:,np) = round([xi;yi]);
end
close(h);
