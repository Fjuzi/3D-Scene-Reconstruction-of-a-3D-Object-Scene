%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Viñolo (marcos.escudero@uam.es)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = createGaussianKernel(nr,nc,sigma)
%% define Gaussian filter
% one or two dimensional signal?
if min(nr,nc) == 1
    dim = max(nr,nc);
else
    dim = min(nr,nc);
end
% define with center (odd size)
if rem(dim,2) == 0
    dim = dim+1;
end
g   = fspecial('gaussian',[dim,1],sigma);
% truncate filter_mask to non-zero elements
g   = g(g~=0);