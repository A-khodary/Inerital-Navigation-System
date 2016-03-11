function [dif] = subtr_ang(ang1,ang2)

temp1 = ang1;
% make theta in [0,2pi)
if(temp1 >0)
    temp1 = temp1 - 2*pi*floor(temp1/(2*pi));
end
if(temp1 <0)
    temp1 = 2*pi - (-temp1 - 2*pi*floor(-temp1/(2*pi)));
end

temp2 = ang2;
% make theta in [0,2pi)
if(temp2 > 0)
    temp2 = temp2 - 2*pi*floor(temp2/(2*pi));
end
if(temp2 < 0)
    temp2 = 2*pi - (-temp2 - 2*pi*floor(-temp2/(2*pi)));
end

dif = temp1 - temp2;
if(dif > pi)
    dif = dif - 2*pi;
end
if(dif < -pi)
    dif = dif + 2*pi;
end