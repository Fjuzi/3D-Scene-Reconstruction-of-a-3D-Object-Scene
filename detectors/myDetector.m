%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Viñolo (marcos.escudero@uam.es)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points,decomposition] =  myDetector(gima,params)
    %% get detections:
    switch params.detector

         case 'LoG_ss'

              SS                  =   doScaleSpaceGivenSigmas(gima,params.sigmas);
              LoG                 = extractLaplacianGivenSigmas(SS,params.sigmas);
              [M,MaxPos,m,MinPos] = MinimaMaxima3D(LoG,1,0,params.npoints,params.npoints);

                %local maxima
              points.Metric          = abs(m);
              points.Location        = [MinPos(:,2),MinPos(:,1)];
              points.SignOfLaplacian = -1;

              points.sid             = MinPos(:,3);
                %localm inima
              points.Metric          = [points.Metric;M];
              points.Location        = [points.Location;[MaxPos(:,2),MaxPos(:,1)]];
              points.SignOfLaplacian = 1;
              points.sid             = [points.sid;MaxPos(:,3)];
              
              points.Scale           = params.sigmas(points.sid);
              
              decomposition          = SS;
              
         case 'SURF'
             
             points                 =  detectSURFFeatures(gima); %matlab implentation of surf
             
             decomposition{1}       = gima;
             
        case 'Kaze'
             
             frame                 =  detectKAZEFeatures(gima,'NumOctaves', params.noctaves, 'NumScaleLevels', params.nscales); %matlab implentation of Kaze
             points.Location        = frame.Location;
             points.Scale           =  frame.Scale;
             points.Orientation     =  frame.Orientation;
             points.Metric          =   frame.Metric;
             points.SignOfLaplacian =        zeros(size(frame.Scale));
             points.Count = frame.Count;
             decomposition{1}       = gima;     
             
        case 'SIFT'
             % see additional detection methods and options provided by vlfeat in help vl_covdet
              [frame]= ...
              vl_sift(single(gima),'Octaves',params.noctaves,'Levels',params.nscales); 
            % parse to common structure
               points.Location        = [frame(1,:)',frame(2,:)'];
               points.Scale           =  frame(3,:)';
               points.Orientation     =  frame(4,:)';
               points.Metric          =   0.5.*ones(size( frame(4,:)'));% unknown
               points.SignOfLaplacian =        zeros(size(frame(4,:)'));     % unknown
               
               decomposition{1}       = gima;    

               
         case 'DoH'
             SS                  =   doScaleSpaceGivenSigmas(gima,params.sigmas);
             DoH                 = extractDeterminantOfHessianGivenSigmas(SS,params.sigmas);
             [M,MaxPos,m,MinPos] = MinimaMaxima3D(DoH,1,0,params.npoints,params.npoints);
             
             Saddle.v = m;Saddle.x = MinPos(:,2);Saddle.y = MinPos(:,1); Saddle.s = MinPos(:,3);
                PoI.v    = M;PoI.x    = MaxPos(:,2);PoI.y = MaxPos(:,1); PoI.s = MaxPos(:,3);

                  %local maxima
              points.Metric          = abs(m);
              points.Location        = [MinPos(:,2),MinPos(:,1)];
              points.SignOfLaplacian = -1;

              points.sid             = MinPos(:,3);
                %localm inima
              points.Metric          = [points.Metric;M];
              points.Location        = [points.Location;[MaxPos(:,2),MaxPos(:,1)]];
              points.SignOfLaplacian = 1;
              points.sid             = [points.sid;MaxPos(:,3)];
              
              points.Scale           = params.sigmas(points.sid);
              
              decomposition          = SS;
                
                
                
                %Without Saddle Points
%               points.Metric          = abs(M);
%               points.SignOfLaplacian = -1;           
%             
%               points.Location        = [PoI.x,PoI.y];
%               points.sid             = MaxPos(:,3);
%               
%               points.Scale           = params.sigmas(points.sid);
%               
%               decomposition          = SS;
             
         case 'anyother' %can be anything KAZE... 
               ...

        otherwise  % SURF

              points                 =  detectSURFFeatures(gima);
             
              decomposition{1}       = gima;

    end
% note that matlab's functions: detectKAZEFeatures, detectSURFFeatures,
% detect<<XXX>>Features,.. allow you to detect using several of the
% handcrafted detectors studied in class.
end