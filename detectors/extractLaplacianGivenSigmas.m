function LoG  = extractLaplacianGivenSigmas(SS,sigmas)

LoG = zeros(size(SS{1},1),size(SS{1},2),numel(SS));
for j=1:numel(SS)
   
    dggxx  = ndgauss([size(SS{j},1) size(SS{j},2)],[sigmas(j),sigmas(j)],'der',[2 0]);
    Lxx    = imfilter(SS{j},dggxx,'symmetric','same');
    
    dggyy  = ndgauss([size(SS{j},1) size(SS{j},2)],[sigmas(j),sigmas(j)],'der',[0 2]);
    Lyy    = imfilter(SS{j},dggyy,'symmetric','same');
    
    normf      = sigmas(j).^(2);% t
    LoG(:,:,j) = normf.*(Lxx + Lyy);
    
end