function [L] = PostScatteringLength(D_L, D_d, O, r, x, a, b)
%Calculates the length through the circle after scattering occured

a_q = (D_d - a)^2 + (x - D_L/2 - b)^2;
b_q = 2*a*(D_d - a) - 2*O*(D_d - a) + 2*b*(x - D_L/2 - b);
c_q = a^2 + b^2 - 2*a*O - r^2 + O^2;
t = (-b_q + sqrt(b_q^2 - 4*a_q*c_q))/(2*a_q); %Quadratic formula
L = 100*sqrt(((D_d - a)*t)^2 + ((x - D_L/2 - b)*t)^2);

end
