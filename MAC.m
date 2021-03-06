function [mew] = MAC(E)
%Returns a Mass Attentuation Coefficient Value
EMeV = E/1000000;

%Interpolation
p = [7.11200E-03  4.076E+02
   8.00000E-03  3.056E+02
   1.00000E-02  1.706E+02 
   1.50000E-02  5.708E+01 
   2.00000E-02  2.568E+01
   3.00000E-02  8.176E+00
   4.00000E-02  3.629E+00 
   5.00000E-02  1.958E+00 
   6.00000E-02  1.205E+00
   8.00000E-02  5.952E-01
   1.00000E-01  3.717E-01
   1.50000E-01  1.964E-01 
   2.00000E-01  1.460E-01 
   3.00000E-01  1.099E-01 
   4.00000E-01  9.400E-02
   5.00000E-01  8.414E-02
   6.00000E-01  7.704E-02 
   8.00000E-01  6.699E-02 
   1.00000E+00  5.995E-02 
   1.25000E+00  5.350E-02 
   1.50000E+00  4.883E-02 
   2.00000E+00  4.265E-02 
   3.00000E+00  3.621E-02 
   4.00000E+00  3.312E-02
   5.00000E+00  3.146E-02
   6.00000E+00  3.057E-02
   8.00000E+00  2.991E-02
   1.00000E+01  2.994E-02
   1.50000E+01  3.092E-02
   2.00000E+01  3.224E-02];

g=zeros(30, 1);
h=zeros(30, 1);
for i=1:30;
    g(i) = p(i, 1);
    h(i) = p(i, 2);
end

mew = interp1(g, h, EMeV, 'spline');

end
