function [L] = PreScatteringLength(CirCen, a, b, r)
%Calculates the length through the circle before scattering occured

t = ((2*CirCen*a) - sqrt((2*CirCen*a)^2 - 4*(a^2 + b^2)*(CirCen^2 - r^2)))/(2*(a^2 + b^2)); %Quadratic formula
L = 100*sqrt((a*(1 - t))^2 + (b*(1 - t))^2); %Attenuation Length [cm]

end
