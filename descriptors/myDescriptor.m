%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Viñolo (marcos.escudero@uam.es)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [features,new_points] =  myDescriptor(points,decomposition,params)
    npoints            = numel(points.Scale);
    points.Orientation = zeros(npoints,1);
    features           = [];
    %% get descriptors:
    switch params.descriptor

         case 'SURF'
             
             % on decomposition ?
             if params.desOnDecom
                 
                 % is decomposition available ?
                 if numel(decomposition) > 1  
             
                     % go through scales/octaves to extract descriptions
                     for j=1:numel(decomposition)

                         points_at_scale   =    points.sid  == j;
                         image_at_scale    =    decomposition{j};

                         if sum(points_at_scale) > 0

                         % cheat extractFeatures by a SURF parser
                         parser                   = SURFPoints(points.Location(points_at_scale,:));
                         % extractFeatures on detected scale
                         [localfeatures, valid_points] = extractFeatures(image_at_scale,parser);

                             % recover orientation
                             if valid_points.Count ~= parser.Count
                                 sprintf('Unable to recover orientation, do not extract points on boundaries with MinimaMaxima3D!');
                             else
                                 points.Orientation(points_at_scale) = valid_points.Orientation;
                             end

                             % store features
                             if isempty(features) % reserve memory when descriptor's length is known
                                 features = rand(npoints,size(localfeatures,2));
                             end
                             features(points_at_scale,:) = localfeatures;% ignore the warning

                         end 

                     end

                     % establish SURF as detected points (for representation)
                     new_points  = ...
                     SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                points.Metric,'Orientation',points.Orientation,...
                                'SignOfLaplacian',points.SignOfLaplacian);
                            
                 else
                     sprintf('Problem: No decomposition is available, changing to standard SURF')
                     
                     switch class(points)
                
                     case {'SURFPoints', 'vision.internal.SURFPoints_cg',}
                         % donothing
                     otherwise    % adapt detected points to SURF descriptor
                         points  = ...
                                   SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                   points.Metric,'SignOfLaplacian',points.SignOfLaplacian);
                     end
             
                     [features, new_points] = extractFeatures(decomposition{1},points);
                 end
                 
             else
                 
                 switch class(points)
                
                     case {'SURFPoints', 'vision.internal.SURFPoints_cg',}
                         % donothing
                     otherwise    % adapt detected points to SURF descriptor
                         points  = ...
                                   SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                   points.Metric,'SignOfLaplacian',points.SignOfLaplacian);
                     
                 end
             
                 [features, new_points] = extractFeatures(decomposition{1},points);
                
             end
             % endof SURF
          case 'SIFT'
             
             % on decomposition ?
             if params.desOnDecom
                 
                 % is decomposition available ?
                 if numel(decomposition) > 1  
             
                     % go through scales/octaves to extract descriptions
                     for j=1:numel(decomposition)

                         points_at_scale   =    points.sid  == j;
                         image_at_scale    =    decomposition{j};

                         if sum(points_at_scale) > 0

                         % cheat vl_covdet by a SIFT parser
                         frame      = zeros(3,sum(points_at_scale));
                         frame(1,:) = points.Location(points_at_scale,1);
                         frame(2,:) = points.Location(points_at_scale,2);
                         frame(3,:) = ones(size(frame(2,:)));
                         
                         [localpoints,localfeatures]  = vl_covdet(single(image_at_scale), 'Frames', frame,...
                                                        'EstimateOrientation', true, 'MaxNumOrientations',1); 
                         % Max number of orientations has been set to 1 to
                         % ease feature allocation (may be improved).
                         
                         % store features
                           localfeatures = localfeatures';
                           if isempty(features) % reserve memory when descriptor's length is known
                                 features = rand(npoints,size(localfeatures,2));
                           end
                           features(points_at_scale,:) = localfeatures;% ignore the warning
                         
                         % recover orientation
                         points.Orientation(points_at_scale) = atan2(localpoints(6,:),localpoints(5,:))';

                         end 

                     end
      
                 else
                     sprintf('Problem: No decomposition is available, changing to standard SIFT')
                     
                     % cheat vl_covdet by a SIFT parser
                     frame      = zeros(3,npoints);
                     frame(1,:) = points.Location(:,1);
                     frame(2,:) = points.Location(:,2);
                     frame(3,:) = points.Scale;
                         
                     [auxpoints,features,addinfo]  = vl_covdet(single(decomposition{1}), 'Frames', frame,...
                                             'EstimateAffineShape', false, 'EstimateOrientation', true, 'MaxNumOrientations',1); 
                    % Max number of orientations has been set to 1 to
                    % ease feature allocation (may be improved).
                    
                    features = features';
                                         
                     
                     % recover score if possible
                     if sum(isnan(addinfo.laplacianScaleScore')) == 0
                     points.Metric = addinfo.laplacianScaleScore';
                     end
                         
                     % recover orientation
                     points.Orientation = atan2(auxpoints(6,:),auxpoints(5,:))';
                     
                 end
                 
             else
                 
                     % cheat vl_covdet by a SIFT parser
                     frame      = zeros(3,npoints);
                     frame(1,:) = points.Location(:,1);
                     frame(2,:) = points.Location(:,2);
                     frame(3,:) = points.Scale;
                         
                     [auxpoints,features,addinfo]  = vl_covdet(single(decomposition{1}), 'Frames', frame,...
                                             'EstimateOrientation', true, 'MaxNumOrientations',1); 
                     % Max number of orientations has been set to 1 to
                     % ease feature allocation (may be improved).
                     
                     features = features';
                                         
                     % recover score if possible
                     if sum(isnan(addinfo.laplacianScaleScore')) == 0
                     points.Metric = addinfo.laplacianScaleScore';
                     end
                         
                     % recover orientation
                     points.Orientation = atan2(auxpoints(6,:),auxpoints(5,:))';
                
             end
             
             % establish SURF as detected points (for representation)
             try
             new_points  = ...
                        SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                        points.Metric,'Orientation',points.Orientation,...
                        'SignOfLaplacian',points.SignOfLaplacian);
             catch
                 keyboard;
             end
             % endofSIFT
             
          case 'DSP-SIFT'
             
             % on decomposition ?
             if params.desOnDecom
                 
                 sprintf('Problem: DSP-SIFT has not been yet coded to operate on decompositions \n Extracting DSP-features on SIFT scale-space')
                 
             end
              
              % cheat vl_covdet by a SIFT parser
               frame      = zeros(3,npoints);
               frame(1,:) = points.Location(:,1);
               frame(2,:) = points.Location(:,2);
               frame(3,:) = points.Scale;
                         
               [localframes,~]  = vl_covdet(single(decomposition{1}), 'Frames', frame,...
                                                        'EstimateOrientation', true, 'MaxNumOrientations',1); 
                                                    
               [localframes, features] = vl_dspcovsift(localframes,decomposition{1},params);
                features = features';
             
                % recover orientation
                points.Orientation = atan2(localframes(6,:),localframes(5,:))';
                 
                % establish SURF as detected points (for representation)
                new_points  = ...
                            SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                            points.Metric,'Orientation',points.Orientation,...
                            'SignOfLaplacian',points.SignOfLaplacian);

                % endofDSP-SIFT

         case 'Kaze'
             
             if numel(decomposition) > 1  
                 for j=1:numel(decomposition)
                     points_at_scale   =    points.sid  == j;
                     image_at_scale    =    decomposition{j};
                     if sum(points_at_scale) > 0
                          parser = KAZEPoints(points.Location(points_at_scale,:));
                         [localfeatures, valid_points] = extractFeatures(image_at_scale,parser);
                     end
                     
                     
                             if valid_points.Count ~= parser.Count
                                 sprintf('Unable to recover orientation, do not extract points on boundaries with MinimaMaxima3D!');
                             else
                                 points.Orientation(points_at_scale) = valid_points.Orientation;
                             end

                             % store features
                             if isempty(features) % reserve memory when descriptor's length is known
                                 features = rand(npoints,size(localfeatures,2));
                             end
                             features(points_at_scale,:) = localfeatures;% ignore the warning
                 end
                 
             else
                parser = KAZEPoints(points.Location, 'Scale',points.Scale, 'Metric',...
                            points.Metric, 'Orientation',points.Orientation);
               [features, new_points] = extractFeatures(decomposition{1},parser);
                 
             end
               
        otherwise  % SURF
            
             % on decomposition ?
             if params.desOnDecom
                 
                 % is decomposition available ?
                 if numel(decomposition) > 1  
             
                     % go through scales/octaves to extract descriptions
                     for j=1:numel(decomposition)

                         points_at_scale   =    points.sid  == j;
                         image_at_scale    =    decomposition{j};

                         if sum(points_at_scale) > 0

                         % cheat extractFeatures by a SURF parser
                         parser                   = SURFPoints(points.Location(points_at_scale,:));
                         % extractFeatures on detected scale
                         [localfeatures, valid_points] = extractFeatures(image_at_scale,parser);

                             % recover orientation
                             if valid_points.Count ~= parser.Count
                                 sprintf('Unable to recover orientation, do not extract points on boundaries with MinimaMaxima3D!');
                             else
                                 points.Orientation(points_at_scale) = valid_points.Orientation;
                             end

                             % store features
                             if isempty(features) % reserve memory when descriptor's length is known
                                 features = rand(npoints,size(localfeatures,2));
                             end
                             features(points_at_scale,:) = localfeatures;% ignore the warning

                         end 

                     end

                     % establish SURF as detected points (for representation)
                     new_points  = ...
                     SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                points.Metric,'Orientation',points.Orientation,...
                                'SignOfLaplacian',points.SignOfLaplacian);
                            
                 else
                     sprintf('Problem: No decomposition is available, changing to standard SURF');
                     
                     switch class(points)
                
                     case {'SURFPoints', 'vision.internal.SURFPoints_cg',}
                         % donothing
                     otherwise    % adapt detected points to SURF descriptor
                         points  = ...
                                   SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                   points.Metric,'SignOfLaplacian',points.SignOfLaplacian);
                     end
             
                     [features, new_points] = extractFeatures(decomposition{1},points);
                 end
                 
             else
                 
                 switch class(points)
                
                     case {'SURFPoints', 'vision.internal.SURFPoints_cg',}
                         % donothing
                     otherwise    % adapt detected points to SURF descriptor
                         points  = ...
                                   SURFPoints(points.Location,'Scale',points.Scale,'Metric',...
                                   points.Metric,'SignOfLaplacian',points.SignOfLaplacian);
                     
                 end
             
                 [features, new_points] = extractFeatures(decomposition{1},points);
                
             end
 

    end
% note that matlab's function: extractFeatures,.. allows you to describe using several of the
% handcrafted descriptors mentioned in class.
end

%% Adapted from
% The DSP-SIFT library
% Copyright (c) 2014-2015, Jingming Dong
% All rights reserved.
function [f_out, d_out] = vl_dspcovsift(frames,ima,params)

    ns     =     params.dsp.ns;
    sc_min = params.dsp.sc_min;
    sc_max = params.dsp.sc_max;
    
    % Sample scales around detections
    nf = size(frames, 2);
    f = zeros(6, nf, ns);
    cnt_sc = 0;
    for sc = linspace(sc_min, sc_max, ns)
        cnt_sc = cnt_sc + 1;
        f([1 2], :, cnt_sc) = frames([1 2], :);
        f(3:6, :, cnt_sc)   = sc * frames(3:6,:);
    end
    
    % Compute un-normalized SIFT at each scales
    f = reshape(f, [6, nf*ns, 1]);
    
    f_size = zeros(1,nf*ns);
    for i=1:nf*ns
        temp = f(:,i);
        E = reshape(temp(3:6),2,2);        
        f_size(i) = det(E);
    end
    
    [~, assign] = sort(f_size);
    f_sort = f(:, assign);
    [~, assign_back] = sort(assign);
    
    [f_out,d]  = vl_covdet(single(ima), 'Frames', f_sort,...
                                        'EstimateAffineShape',false,'EstimateOrientation', false, 'MaxNumOrientations',1); 
                                         
    f_out = f_out(:, assign_back);
    d = d(:, assign_back);
    
    % Aggregate and nomarlize
    f_out = reshape(f_out, [6, nf, ns]);
    f_out = f_out(:, :, 1);
    f_out(3:6,:) = frames(3:6,:);
    d = reshape(d, [128, nf, ns]);
    
    d_out = mean(single(d), 3);

    % normalization...
    d_out = d_out./repmat(sqrt(sum(d_out.^2)),[size(d_out,1),1]);
    
end