%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Viñolo (marcos.escudero@uam.es)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SS=doScaleSpaceGivenSigmas(I,sigmas)

%% Initialize scale-space structure
[nr,nc]     = size(I);
nscales     = numel(sigmas);
SS          = cell(1, nscales);
SS{1}       = double(I);
%% obtain scale-space representation guided by the Gaussian filter
sigmas = [0, sigmas];
for j = 2:1:(nscales)
    g     = createGaussianKernel(nr,nc,sigmas(j));
    SS{j} =  imfilter(SS{1},g,'symmetric','same');
end

end
