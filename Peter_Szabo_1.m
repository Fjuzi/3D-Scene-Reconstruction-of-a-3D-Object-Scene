close all;
clear all;

%number of points
Np = 9;
%width of the image
widhtP = 360;
% 0: no display
display = 0;
[coords, ima_pattern]= get_real_points_checkerboard_vmmc(Np, widhtP, display);
%load('xy_origin.mat');

load('marked_targets_BIG.mat');

%% BIG images - read the images (you can speicfy yourself, needed at lest 6-7 images to work)
img1 = imread('20200506_133024.jpg');
% xy_target = get_user_points_vmmc(img1);
% marked_targets(:,:,1) = xy_target;
% 
img2 = imread('20200506_133031.jpg');
% xy_target = get_user_points_vmmc(img2);
% marked_targets(:,:,2) = xy_target;
% 
img3 = imread('20200506_133039.jpg');
% xy_target = get_user_points_vmmc(img3);
% marked_targets(:,:,3) = xy_target;
% 
img4 = imread('20200506_133047.jpg');
% xy_target = get_user_points_vmmc(img4);
% marked_targets(:,:,4) = xy_target;
% 
img5 = imread('20200506_133027.jpg');
% xy_target = get_user_points_vmmc(img5);
% marked_targets(:,:,5) = xy_target;
% 
img6 = imread('20200506_133035.jpg');
% xy_target = get_user_points_vmmc(img6);
% marked_targets(:,:,6) = xy_target;
% 
img7 = imread('20200506_133044.jpg');
% xy_target = get_user_points_vmmc(img7);
% marked_targets(:,:,7) = xy_target;
%save('marked_targets_BIG.mat','marked_targets')
%% Convert to the required format
img{1} = img1;
img{2} = img2;
img{3} = img3;
img{4} = img4;
img{5} = img5;
img{6} = img6;
img{7} = img7;
montage(img);
%% Create homography
for j = 1:7
    v = homography_solve_vmmc(coords', marked_targets(:,:,j));
    [H_refined, rperr] = homography_refine_vmmc(coords', marked_targets(:,:,j), v);
    Hex3{j} = H_refined;
    T = maketform('projective', H_refined');
    tr_ima = imtransform(ima_pattern,T,'XData',[1 size(img1,2)], 'YData',[1 size(img1,1)]);
    %transformed_images(:,:,1) = tr_ima;
    %imshow(tr_ima);
    %ause;
end
% Caculate the parameters of the inner matrix
h = cell(1,3);
h{1} = Hex3{1};
h{2} = Hex3{2};
h{3} = Hex3{3};
A1 = internal_parameters_solve_vmmc(h);
mcos3 = radtodeg(acos (-1 * A1(1,2)/A1(1,1)));
[R3, T3] = external_parameters_solve_vmmc(A1, h);


h{4} = Hex3{4};
A2 = internal_parameters_solve_vmmc(h);
mcos4 = radtodeg(acos (-1 * A2(1,2)/A2(1,1)));
[R4, T4] = external_parameters_solve_vmmc(A2, h);


h{5} = Hex3{5};
A3 = internal_parameters_solve_vmmc(h);
mcos5 = radtodeg(acos (-1 * A3(1,2)/A3(1,1)));
[R5, T5] = external_parameters_solve_vmmc(A3, h);


h{6} = Hex3{6};
A4 = internal_parameters_solve_vmmc(h);
mcos6 = radtodeg(acos (-1 * A4(1,2)/A4(1,1)));
[R6, T6] = external_parameters_solve_vmmc(A4, h);

h{7} = Hex3{7};
A5 = internal_parameters_solve_vmmc(h);
mcos7 = radtodeg(acos (-1 * A5(1,2)/A5(1,1)));
[R7, T7] = external_parameters_solve_vmmc(A5, h);
%alpha = radtodeg(atan( A(1,2) / A(1,1)));

%save('marked_targets_SMALL.mat','marked_targets')


%% Small images: (this part is only for the experiment for the smaller point, it wont be used later on
close all;
clear all;

Np = 9;
widhtP = 290;
display = 0;
[coords, ima_pattern]= get_real_points_checkerboard_vmmc(Np, widhtP, display);
%load('xy_origin.mat');

load('marked_targets_SMALL.mat');


img1 = imread('20200506_150540.jpg');
%xy_target = get_user_points_vmmc(img1);
%marked_targets(:,:,1) = xy_target;

img2 = imread('20200506_150553.jpg');
% xy_target = get_user_points_vmmc(img2);
% marked_targets(:,:,2) = xy_target;

img3 = imread('20200506_150604.jpg ');
% xy_target = get_user_points_vmmc(img3);
% marked_targets(:,:,3) = xy_target;

img4 = imread('20200506_150534.jpg');
% xy_target = get_user_points_vmmc(img4);
% marked_targets(:,:,4) = xy_target;

img5 = imread('20200506_150545.jpg');
% xy_target = get_user_points_vmmc(img5);
% marked_targets(:,:,5) = xy_target;

img6 = imread('20200506_150558.jpg');
% xy_target = get_user_points_vmmc(img6);
% marked_targets(:,:,6) = xy_target;




