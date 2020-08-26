function [coords ima_pattern]= get_real_points_checkerboard_vmmc(Np, widhtP, display)
% GET_REAL_POINTS_CHECKERBOARD_VMMC(Np, widhtP, display)  returns a list of Np points 
% automatically selected for a checkerboard pattern with real dimensions 
% defined by the witdh of the checkerboard pattern (widhtP)
%
% Use display equal to 1 to show the order for the returned point coordinates.
%
% In practical situations, 9 points are required for convergence due to
% numerical reasons. Moreover, these points have to be distributed
% throughout the image as much as possible.
%
%
% Returns a 2xNp matrix representing the point coordinates (only 2D XY coordinates
%   are returned as Z=0 is assumed for the pattern) and the synthetic image pattern. 
%
% Juan C. SanMiguel, Universidad Autonoma de Madrid, November 2012

%compute coordinates of pattern points
step = (sqrt(Np)-1)/2;
m = 1;
for j = widhtP : -widhtP/(2*step) : 0
    for i = 0 : widhtP/(2*step) : widhtP
        %coords(m, :) = [i, j, 0, 1]; %3D coords with Z = 0
        coords(m, :) = [i, j];%2D coords
        m = m + 1;
    end
end

ima_pattern=checkerboard(round(widhtP/8)) > 0.5;

if display == 1
    %plot points over the pattern and its order
    x = coords(:,1); x(x==0)=1;
    y = coords(:,2); y=abs(y);y(y==0)=1;
        
    figure; imshow(ima_pattern);
    hold on; plot(x,y,'r*');
    for j = 1:numel(x)
        text(x(j),y(j),sprintf('  %d',j), 'Color', [1 0 0])
    end
end