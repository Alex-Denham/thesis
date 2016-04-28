%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Input Variables  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
Res = 10;         %Resolution of model
r = 0.005;        %Radius of sphere
E1 = 225E3;       %Emitter energy [eV]
D_L = 0.04;       %Length of detector [m]
D_d = 0.1;        %Distance between emitter and detector [m]
O = 0.09;         %Distance from source to the centre of the circle [m]
I_0 = 1E12;       %Source Intensity [Photons]
No_PP = 200;      %Number of plot points

%%%%%%%%%%%%%%
%%%  Code  %%%
%%%%%%%%%%%%%%

%Constants
e = 1.602176620898E-19;     %Elementary charge
Na = 6.02214086E23;         %Avogadro's
E1_J = E1*e;                %Conversion to Joules [J]
MAC_E1 = MAC(E1);           %Mass Attenuation Coefficent for given energy and material
rho = 7.874;                %Density of Iron [g/cm^3]
ES = (2*r)/Res;             %Size of element [m]
PS = D_L/No_PP;             %Pixel size
AM = 55.845;                %Atomic mass [Da]
N = (1000000*rho*ES*Na)/AM; %Number of atoms in a single element
phi = 2*atan(D_L/(2*D_d));  %Total angle of photon distribution from source.

%Initiation
x = linspace(0, D_L, No_PP);          %Location on detector (j direction) [m]
[yR, yC, yN] = deal(zeros(1, No_PP)); %Better for memory allocation and speed

%scattering
for i = 1:Res
    a = i*ES - ES/2 + O - r; %location of element centre i direction [m]
    for j = 1:Res        
        b = j*ES - ES/2 - r; %location of element centre j direction [m]        
        zeta_1 = abs(atan((abs(b)+ES/2)/(a-ES/2))-atan((abs(b)-ES/2)/(a-ES/2)))/phi; % Percentage light intensity that hits voxel (j direction)          
        AL1 = PreScatteringLength(O, a, b, r);  %Attenuation length before scattering [cm]
        if sqrt(b^2 + (a-O)^2) < r
            for k = 1:Res
                c = k*ES - ES/2 - r;
                xi_1 = atan(abs(c)/a);       %Vertical angle incidence (for attenuation)
                xi_2 = atan(abs(c)/(D_d-a)); %Vertical angle scattered (for attenuation)
                zeta_2 = abs(atan((abs(c)+ES/2)/(a-ES/2))-atan((abs(c)-ES/2)/(a-ES/2)))/phi; % Percentage light intensity that hits voxel (vertical)
                I_1 = I_0*zeta_1*zeta_2;     %Incident on the element
                psi = abs(atan((abs(c)+PS/2)/(D_d-a))-atan((abs(c)-PS/2)/(D_d-a)));
                for n = 1:No_PP                
                    u = [a, b, c];
                    v = [D_d-a, x(n)-D_L/2 - b, -c];
                    theta= acos(dot(u, v)/(norm(u)*norm(v))); %Scatter angle
                    AL2 = PostScatteringLength(D_L, D_d, O, r, x(n), a, b); %Attenuation Length from element to detector [cm]                
                    gamma = abs(atan((abs(x(n)-b-D_L/2)+(PS/2))/(D_d-a))-atan((abs(x(n)-b-D_L/2)-(PS/2))/(D_d-a))); %angle of scatter for a given pixel
                    I_2 = I_1 * exp(-MAC_E1*rho*((AL1/cos(xi_1))+(AL2/cos(xi_2)))); %Intensity after attenuation
                    yR(n) = yR(n) + N*I_2*RayleighFF(theta, E1_J)*(gamma*psi);
                    yC(n) = yC(n) + N*I_2*ComptonISF(theta, E1_J)*(gamma*psi);
                end
            end
        end
    end
    Progress = i*100/Res
end

%Direct
zeta_3 = (2*atan(PS/(2*D_d)))/phi; %Ratio of light that hits pixel (depth)
for n = 1:No_PP
    I_3 = I_0*zeta_3*abs(atan((abs(x(n)-(D_L/2))+PS)/D_d)-atan((abs(x(n)-(D_L/2))-PS)/D_d))/phi; %Intensity on each pixel
    AL3 = LengthThroughCircle(D_L, D_d, O, r, x(n)); %Non scattered attenuation length
    if isreal(AL3) == 1
        yN(n) = I_3 * exp(-MAC_E1*rho*AL3); %Insensity at detector
    else
        yN(n) = I_3;
    end
end

yT = yR + yC + yN; %Total incident on the detector pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Detector image from a  %%%
%%%%%  single projection  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PP = 1:No_PP;

figure
plot(PP, yR, 'Color', [0 0 0]);
title('Rayleigh Scattering Line Plot')
xlabel('Pixel Number')
ylabel('Intensity [Photons]')

figure
plot(PP, yC, 'Color', [0 0 0]);
title('Compton Scattering Line Plot')
xlabel('Pixel Number')
ylabel('Intensity [Photons]')

figure
QT1 = plot(PP, yN, 'Color', [0 0 0]);
hold on
QT2 = plot(PP, yT, 'Color', [0.6 0.6 0.6]);
hold off
title('Total Intensity')
xlabel('Pixel Number')
ylabel('Intensity [Photons]')
legend([QT1, QT2], {'Direct', 'Direct + Scattered'})

%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Reconstruction  %%%
%%%%%%%%%%%%%%%%%%%%%%%%

d = round(D_d/PS);
No_rot = 360;
[fR, fC, fN, fT] = deal(zeros(No_PP, No_rot));

for z = 1:No_rot
    fR(:, z) = yR;
    fC(:, z) = yC;
    fN(:, z) = yN;
    fT(:, z) = yT;
end

iN = ifanbeam(fN, d, 'FanSensorGeometry', 'line');
iN = iN/(min(min(iN))); %Scaling

iT = ifanbeam(fT, d, 'FanSensorGeometry', 'line');
iT = iT/(min(min(iT))); %Scaling

clims = [-1 0.5]; %Sets limits for colour range
IS = size(iN);    %Size of reconstructed image

figure
imagesc(-iN, clims)
colormap(gray)
axis equal
axis([0 IS(1) 0 IS(2)])
xlabel('Pixel Number')
ylabel('Pixel Number')
title('Reconstruction Direct Only')

figure
imagesc(-iT, clims)
colormap(gray)
axis equal
axis([0 IS(1) 0 IS(2)])
xlabel('Pixel Number')
ylabel('Pixel Number')
title('Reconstruction Direct and Scattered')

idif = 100*(iT - iN); %A plot of the difference between data sets
figure
imagesc(idif)
axis equal
title('Difference Comparison')

%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Edge detection  %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Initially set all values at 1
ED3 = ones(IS(1), IS(2));
ED4 = ones(IS(1), IS(2));

%Loop that determines whether a pixel is on the edge of the circle.
for i = 2:(IS(1)-1)
    for j = 2:(IS(2)-1)
        if (iN((i-1),(j-1))<0 || iN((i+1),(j-1))<0 || iN((i-1),(j+1))<0 || iN((i+1),(j+1))<0 ...
                 || iN((i),(j-1))<0 || iN((i),(j+1))<0 || iN((i-1),(j))<0 || iN((i+1),(j))<0) ...
                 && (iN((i-1),(j-1))>0 || iN((i+1),(j-1))>0 || iN((i-1),(j+1))>0 || iN((i+1),(j+1))>0 ...
                 || iN((i),(j-1))>0 || iN((i),(j+1))>0 || iN((i-1),(j))>0 || iN((i+1),(j))>0)
            ED3(i, j) = 0;
        end
        if (iT((i-1),(j-1))<0 || iT((i+1),(j-1))<0 || iT((i-1),(j+1))<0 || iT((i+1),(j+1))<0 ...
                 || iT((i),(j-1))<0 || iT((i),(j+1))<0 || iT((i-1),(j))<0 || iT((i+1),(j))<0) ...
                 && (iT((i-1),(j-1))>0 || iT((i+1),(j-1))>0 || iT((i-1),(j+1))>0 || iT((i+1),(j+1))>0 ...
                 || iT((i),(j-1))>0 || iT((i),(j+1))>0 || iT((i-1),(j))>0 || iT((i+1),(j))>0)
            ED4(i, j) = 0;
        end
    end
end

%Plots of the edge
figure, imshow(ED3), axis equal, title('Edge Detection (Direct)')
figure, imshow(ED4), axis equal, title('Edge Detection (Direct+Scattered)')

%Turning the edge pixels into data points
mN = 1;
mT = 1;
xfN=0;
yfN=0;
xfT=0;
yfT=0;
for i = 1:IS(1)
    for j = 1:IS(2)
        if ED3(i, j) == 0
            xfN(mN)=j;
            yfN(mN)=i;
            mN=mN+1;
        end
        if ED4(i, j) == 0
            xfT(mT)=j;
            yfT(mT)=i;
            mT=mT+1;
        end
    end
end

%Using the edge data points, a circle of best fit is calculated
[xcN, ycN, RN] = circfit(xfN,yfN);
[xcT, ycT, RT] = circfit(xfT,yfT);

th = linspace(0,2*pi,100);
xeN = RN*cos(th)+xcN; yeN = RN*sin(th)+ycN;
xeT = RT*cos(th)+xcT; yeT = RT*sin(th)+ycT;

figure
plot(xfN,yfN,'b.', 'Color', [0.6 0.6 0.6])
hold on
plot(xeN,yeN, 'Color', [0 0 0])
hold off
axis equal
axis([0 IS(1) 0 IS(2)])
xlabel('Pixel Number')
ylabel('Pixel Number')
title('Circle Best Fit')

%Final calculations for the radii
RadiusN = PS*RN*(O/D_d);
RadiusT = PS*RT*(O/D_d);
