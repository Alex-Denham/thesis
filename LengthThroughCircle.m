function [L] = LengthThroughCircle(D_L, D_d, O, r, x)
t1 = ((2*O*D_d) + sqrt((2*O*D_d)^2 - 4*(D_d^2 + (x-(D_L/2))^2)*(O^2 - r^2)))/(2*(D_d^2 + (x-(D_L/2))^2));
t2 = ((2*O*D_d) - sqrt((2*O*D_d)^2 - 4*(D_d^2 + (x-(D_L/2))^2)*(O^2 - r^2)))/(2*(D_d^2 + (x-(D_L/2))^2));

S1 = sqrt((D_d - D_d*t1)^2 + ((x-(D_L/2)) - (x-(D_L/2))*t1)^2);
S2 = sqrt((D_d - D_d*t2)^2 + ((x-(D_L/2)) - (x-(D_L/2))*t2)^2);

L = 100*(S2 - S1); %Non scattered attenuation length

end
