close all;
clear all;
clc;

MaxRatio     =   0.6;
Metric       =  'SSD';
SCALE = 0.2;

params.Directory    = fullfile('../Data/Part23');
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);
ima{numImages}           = [];

for j = 1:numImages
ima{j}      =       imresize(readimage(params.Scene, j),SCALE);
end


load('points.mat');
load('features.mat');
point_matrix = n_view_matching(points, features, ima, MaxRatio, Metric);
%%
% point_matrix = q_r;
% [h, npoints, c] = size(point_matrix);
% q_r = q_r;
% numImages = c;
% q_2cams = point_matrix(:,:,[1,8]);
% [F, P_2cam,Q_2cam,q_2cam] = MatFunProjectiveCalib(q_2cams);
%%
[h, npoints, c] = size(point_matrix);
q_r = ones(h+1,npoints,c);
q_r(1:2,:,:) = point_matrix;
q_2cams = ones(h+1,npoints,2);
q_2cams(1:2,:,:) = point_matrix(:,:,[1,7]);
[F, P_2cam,Q_2cam,q_2cam] = MatFunProjectiveCalib(q_2cams);

['Mean reprojection error from 2 cameras   = ' num2str( ErrorRetroproy(q_2cams,P_2cam,Q_2cam)/2 )]
draw_reproj_error(q_2cams,P_2cam,Q_2cam);
%% show matched images:
visualize1 = permute(point_matrix(:,:,1),[2, 1]);
visualize2 = permute(point_matrix(:,:,7),[2, 1]);
figure(numImages+1); showMatchedFeatures(ima{1},ima{7},visualize1,visualize2);
legend('matched points 1','matched points 7');
title('Matched point between image 1 and image 7');
%% resectioning:
P = zeros(3,4,numImages);
for i = 1:numImages
    [Pi,cost]= PDLT_NA(q_r(:,:,i),Q_2cam);
    P(:,:,i) = Pi; 
end

['Resectioning   = ' num2str( ErrorRetroproy(q_r,P,Q_2cam)/2 )]
draw_reproj_error(q_r,P,Q_2cam);

%% Bundle adjustment
vp = ones(npoints,numImages);
[P,X3d,xc] = BAProjectiveCalib(q_r,P,Q_2cam,vp);
%%
['Bundle adjusted   = ' num2str( ErrorRetroproy(q_r,P,X3d)/2 )]
draw_reproj_error(q_r,P,X3d);

%%
P_final = {P(:,:,1),P(:,:,7)};
F_corrected = vgg_F_from_P(P_final);
%%
load('A3.mat');
%A3 = A3 * 0.2;
A3 = A3 * SCALE;
A3(3,3) = 1;
%A3 = K(:,:,1);
E = A3' * F_corrected * A3;
[R_est,T_est] = factorize_E(E);

Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);

Rcam(:,:,1,1) = eye(3);
Rcam(:,:,2,1) = R_est(:,:,1); 
Tcam(:,2,1) = T_est;

Rcam(:,:,1,2) = eye(3);
Rcam(:,:,2,2) = R_est(:,:,1); 
Tcam(:,2,2) = -1 * T_est;

Rcam(:,:,1,3) = eye(3);
Rcam(:,:,2,3) = R_est(:,:,2); 
Tcam(:,2,3) = T_est;

Rcam(:,:,1,4) = eye(3);
Rcam(:,:,2,4) = R_est(:,:,2); 
Tcam(:,2,4) = -1 * T_est;
%%
close all;
npoints = size(q_r,2);
Q_euc = zeros(4,npoints,2); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;
%K = A3;
K = ones(3,3,2) .* A3;
q_est = xc;
Qeucs = zeros(4,npoints,4);
reprojectedPoints = zeros(3,npoints,2,4);
for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc = TriangEuc(squeeze(Rcam(:,:,2,sol)),squeeze(Tcam(:,2,sol)),K,q_est);
    Qeucs(:,:,sol) = Q_euc; 
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc, K, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
    % Compute the projection matrices from K, Rcam, Tcam
    for k=1:2
       P_euc(:,:,k) = [squeeze(Rcam(:,:,k,sol)), squeeze(Tcam(:,k,sol))];
    end
    
    % Obtain the re-projected points q_rep
    q_rep = P_euc(:,:,1) * Q_euc;
    q_rep(:,:,2) = P_euc(:,:,2) * Q_euc;
    reprojectedPoints(:,:,:,sol) = q_rep;
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
      daspect([1, 1, 1]);
      pbaspect([1, 1, 1]);
      axis([-1000, 1000, -1000, 1000]);
    end
end
%%
['Residual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(reprojectedPoints(:,:,:,1),P_euc,Qeucs(:,:,1))/2 )]
draw_reproj_error(q_r,P,X3d);
%draw_3Dcube(Q_euc(:,:,1));
%%
Q = Qeucs(:,:,1);
%
x = Q(1,:);
y = Q(2,:);
z = Q(3,:);
%scatter3(x,y,z);
scatter(x,y);
X = x;
X(2,:) = y;
%surf(x,y,z);

[idx,C, sum, D] = kmeans(X',2);
%idx((D(:,1) + D(:,2)) > 0.55) = 0;
idx((D(:,1) + D(:,2)) > 0.7) = 0;
% figure;
% plot(X(1,idx==1),X(2,idx==1),'r.','MarkerSize',12)
% %egyik = [X(1,idx==1), X(2,idx==1), z(idx==1)];
% hold on
% plot(X(1,idx==2),X(2,idx==2),'b.','MarkerSize',12)
% %masik = [X(1,idx==2), X(2,idx==2), z(idx==2)];
%
x1 = X(1,idx==1);
y1 = X(2,idx==1);
z1 = z(idx==1);

x2 = X(1,idx==2);
y2 = X(2,idx==2);
z2 = z(idx==2);


stem3(x,y,z)
grid on
xv = linspace(min(x1), max(x1), 20);
yv = linspace(min(y1), max(y1), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x1,y1,z1,X,Y);
figure(2)
hold on;
surf(X, Y, Z);

xv2 = linspace(min(x2), max(x2), 20);
yv2 = linspace(min(y2), max(y2), 20);
[X2,Y2] = meshgrid(xv2, yv2);
Z2 = griddata(x2,y2,z2,X2,Y2);
surf(X2, Y2, Z2);

draw_scene(Q, K, Rcam(:,:,:,1), Tcam(:,:,1));
%%
figure;
hold on;
scatter3(x1,y1,z1, 'r.');
scatter3(x2,y2,z2, 'b.');
draw_scene(Q, K, Rcam(:,:,:,1), Tcam(:,:,1));

% plot(X(1,idx==1),X(2,idx==1),'r.','MarkerSize',12)
% %egyik = [X(1,idx==1), X(2,idx==1), z(idx==1)];
% hold on
% plot(X(1,idx==2),X(2,idx==2),'b.','MarkerSize',12)
%masik = [X(1,idx==2), X(2,idx==2), z(idx==2)];