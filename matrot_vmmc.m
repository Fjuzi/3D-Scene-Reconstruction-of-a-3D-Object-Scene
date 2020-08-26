function [r]=matrot_vmmc(x,unit)

% a) The rotation matrix is calculated from their rotation angles
% b) The rotation angles are calculated from their matrix rotation

% function [r]=matrot(x)
%	unit=0;
%
% function [r]=matrot(x,unit)
% 	unit--> 
%			0: angles in radians (default)
%			1: angles in degrees
% if 'x' is an array then:
% 	x --> angles
% 	r --> rotation matrix
% if 'x' is not an array then:
% 	r --> angles
% 	x --> rotation matrix

% ang(1) --> alfa: rotation round X axis
% ang(2) --> beta: rotation round Y axis
% ang(3) --> gamma: rotation round Z axis
% Positive sense is in the sense of the needles of the clock
% r=Rz*Ry*Rx

% matrot    
%  by Josep Isern González
%  IUSIANI - University of Las Palmas de Gran Canaria
%  January 2002


if (size(x,1)==1) | (size(x,2)==1)
    % ROTATION MATIX CALCULATION
    if nargin==2
        if unit==1
            x=x*pi/180;
        end
    end
    
    sa=sin(x(1));
    ca=cos(x(1));
    sb=sin(x(2));
    cb=cos(x(2));
    sg=sin(x(3));
    cg=cos(x(3));
    
    r=[cb*cg sa*sb*cg-ca*sg ca*sb*cg+sa*sg 0
        cb*sg sa*sb*sg+ca*cg ca*sb*sg-sa*cg 0
        -sb   sa*cb          ca*cb          0
        0     0              0              1];
    
else
   % ROTATION ANGLES CALCULATION
   x2 = x/norm(x);
   g = atan2(x2(2,1),x2(1,1));
   cg = cos(g);
   sg = sin(g);
   b = -atan2(x2(3,1),x2(1,1)*cg+x2(2,1)*sg);
   a = atan2(x2(3,2),x2(3,3));
   if abs(sin(a)*sin(b)*sg+cos(a)*cg-x2(2,2))>1
      x2 = -x/norm(x);
      g = atan2(x2(2,1),x2(1,1));
      cg = cos(g);
      sg = sin(g);
      b = -atan2(x2(3,1),x2(1,1)*cg+x2(2,1)*sg);
      a = atan2(x2(3,2),x2(3,3));
   end
   r = [a,b,g];
   
   if nargin==2
      if unit==1
         r=r*180/pi;
      end
   end
   
end      
