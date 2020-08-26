clear all
close all
clc


params.Directory    = fullfile('../Data/Part23/17');

SCALE = 0.2;
szoveg = '1. - 7. Kaze';
%%Detector parameters
params.detector     =  'Kaze'; %'SURF', 'SIFT', DoH, 'Kaze'
params.nscales      =        10; %5
params.noctaves     =        8; %3
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      500; %300 volt eredetileg

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'Kaze'; % 'SIFT', 'SURF', 'Kaze', 'DSP-SIFT'
params.desOnDecom   =    true; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false; % set to true to avoid orientation estimation. %IT can be avoioded, but harms the dataset

% for DSP-SIFT
params.dsp.ns       =      10;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       3;% largest scale (relative to detection);

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.6; %alpha value for SIFT 0.6 volt eredetileg
params.Metric       =  'SSD'; %explore which metrics are here
% addpaths
addpath(genpath('./detectors/'));
addpath(genpath('./descriptors/'));
addpath(genpath('./toolbox/'));

% preload dataset
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);
ima{numImages}           = [];
points{numImages}        = [];
decomposition{numImages} = [];
features{numImages}      = [];
% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end
% detect & describe
for j = 1:numImages
% Load and convert images %%
%ima{j}      =       readimage(params.Scene, j);
ima{j}      =       imresize(readimage(params.Scene, j),SCALE);
gima        =      im2double(rgb2gray(ima{j})); %convert image grayscale

% PoI Detection %%
%sprintf('Detecting for image: %d',j)
[points{j},decomposition{j}] =  myDetector(gima,params); %inpuit image
%utputs the scale scale representation
% PoI Description %%
%sprintf('Describing for image: %d',j)
[features{j},points{j}]      =  myDescriptor(points{j},decomposition{j},params);

% show detections
figure(j)
imshow(ima{j}); hold on;
plot(points{j},'showOrientation',true);

end

%montage(ima)

% PoI Matching (assumes two images, i.e. numImages == 2) %%
indexPairs       = matchFeatures(features{1},features{2},'MaxRatio',params.MaxRatio,'Metric',params.Metric) ;
matchedPoints{1} = points{1}(indexPairs(:,1));
matchedPoints{2} = points{2}(indexPairs(:,2));
figure(numImages+1); showMatchedFeatures(ima{1},ima{2},matchedPoints{1},matchedPoints{2});
legend('matched points 1','matched points 2');
title(szoveg);

[F_est, varg] = estimateFundamentalMatrix(matchedPoints{1},matchedPoints{2});
[F_est_ransac, varg_ransac] = estimateFundamentalMatrix(matchedPoints{1},matchedPoints{2}, 'Method', 'RANSAC');
%[size(varg,1), sum(varg), sum(varg)/size(varg,1)]
%[size(varg_ransac,1), sum(varg_ransac), sum(varg_ransac)/size(varg_ransac,1)]
%F_est
%% Experiment with the epipoles
fig = vgg_gui_F(ima{1},ima{2},F_est');
%% Visualize points used to create fundamental matrix
matchedUPoints{1} = matchedPoints{1}(varg);
matchedUPoints{2} = matchedPoints{2}(varg);
figure; showMatchedFeatures(ima{1},ima{2},matchedUPoints{1},matchedUPoints{2});
legend('matched points 1','matched points 2');
title(szoveg);

%% Visualize RANSAC points
matchedRPoints{1} = matchedPoints{1}(varg_ransac);
matchedRPoints{2} = matchedPoints{2}(varg_ransac);
figure; showMatchedFeatures(ima{1},ima{2},matchedRPoints{1},matchedRPoints{2});
legend('matched points 1','matched points 2');
title(szoveg);


