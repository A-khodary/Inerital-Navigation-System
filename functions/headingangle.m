function [angle] = headingangle(m)

angle = 0;
% if m(1) < 0
%     angle = pi() - atan(m(2)/m(3));
% end
% if m(1) > 0 & m(2) > 0
%     angle = 2*pi() - atan(m(2)/m(3));
% end
% if m(1) > 0 & m(2) < 0
%     angle =  - atan(m(2)/m(3));
% end
% if m(1) == 0 & m(2) < 0
%     angle =  pi()/2;
% end
% if m(1) == 0 & m(2) > 0
%     angle =  3*pi()/2;
% end
% Para os dados   x_
                   %|\
                   %z y
% if m(2) > 0 
%     angle = pi/2 - atan(m(1)/m(2));
% elseif m(2)<0
%     angle = pi + pi/2 - atan(m(1)/m(2));
% elseif m(2)==0 & m(1)<0
%     angle = pi;
% elseif m(2)==0 & m(1)>0
%     angle = 0;
% end
                  %z y
% Para os dados  x_|/
m(2) = - m(2);
if m(2) > 0 
    angle = pi/2 - atan(m(1)/m(2));
elseif m(2)<0
    angle = pi + pi/2 - atan(m(1)/m(2));
elseif m(2)==0 & m(1)<0
    angle = pi;
elseif m(2)==0 & m(1)>0
    angle = 0;
end

if angle > pi
    angle = -2*pi + angle;
end


