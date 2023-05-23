% The code was provided by André Liemert (https://www.researchgate.net/profile/Andre-Liemert)
% Institut für Lasertechnologien in der Medizin und Meßtechnik an der Universität Ulm, Helmholtzstr. 12, D-89081 Ulm, Germany
% The original file is: [t,R] = DE_time(mua1,mua2,mua3,mus1,mus2,mus3,n1,n2,n3,L1,L2,L3)
% The current file was slightly modified by Aleh Sudakou: added functions "calc_Reff" and "R_fresnel", and the option to convolve with IRF
% This is one of the codes shared on  https://github.com/asudakou/Analyzing_TD-NIRS
% Last updated: 20 May 2023


%% Description
% Generate DTOF for specificed optical properties in a 3 layered medium.
% Returns 'time' in units of ns and the corresponding DTOF in variable 'R'
% 
% PUBLICATIONS:
% A. Liemert, and A. Kienle, "Light diffusion in N-layered turbid media: frequency and time domains," Journal of biomedical optics 15, 025002 (2010).
% DOI:  https://doi.org/10.1117/1.3368682
% 
% A. Liemert, and A. Kienle, "Application of the Laplace transform in time-domain optical spectroscopy and imaging," Journal of biomedical optics 20, 110502 (2015).
% DOI:  https://doi.org/10.1117/1.JBO.20.11.110502
%
% The extrapolated boundary conditions are implemented using equations in publication:
% "Boundary conditions for the diffusion equation in radiative transfer"
% R. C. Haskell, L. O. Svaasand, T. T. Tsay, T. C. Feng, M. McAdams, and B. J. Tromberg,
% Optical Society of America, A 11, 2727–2741 (1994).
% DOI:  https://doi.org/10.1364/JOSAA.11.002727


function [R, time_ns] = DTOF_generate_Liemert(OptProp_all, n, rho, time_s, varargin)


do_this = 0; % This must always be 0
if do_this == 1
    %% Run this section to see Example (Ctrl + Enter)
    OptProp_all_example = [0.01 -1 -1 1 -1 -1 0 0]; % Homogeneous optical properties, Mua=0.01 and Musp=1
    [R_example, time_ns_example] = DTOF_generate_Liemert(OptProp_all_example, 1.4, 30, -1); % Call the current function
    R_example(R_example<0) = 0;
    figure(5); clf
%     semilogy(time_ns, R) % Simple way of plotting DTOF
    [~, ~] = DTOF_filter( R_example, time_ns_example, 0, 0, 0, 0, {5,'red'}); % Better way of plotting DTOF
    clear OptProp_all_example R_example time_ns_example;
end


% INPUTS:
% * OptProp_all :  8 parameters  [Mua1 Mua2 Mua3 Musp1 Musp2 Musp3 L1 L2]    [Units mm-1] and [mm for L1 and L2]
%
% * n           :  Refractive index, e.g. [1.33]  or  [1.33 1.33 1.33]
%
% * rho         :  Source-detector distance in mm, e.g. 30
% 
% * time        :  A vector of time channels. If set to '-1', then will use default 0.1 to 7 ns
% 
% * varargin    :  irf_shifted (Optional). The generated DTOF will be convolved with this IRF
% irf_shifted   :  Photon counts of IRF, which must correspond to the same time channels as "time_s", otherwise need another code for convolution


mua1 = OptProp_all(1); % Absorption and Scattering of 3 layers  [mm-1]
mua2 = OptProp_all(2);
mua3 = OptProp_all(3);
mus1 = OptProp_all(4);
mus2 = OptProp_all(5);
mus3 = OptProp_all(6);
L1   = OptProp_all(7); % Thickness of the first layer
L2   = OptProp_all(8); % Thickness of the second layer
L3 = 60; % Set the third layer thickness to semi-infinite

    if sum(abs([L1, L2])) == 0 % if all L are 0, then it is a homogeneous medium
        L1 = 30;
        L2 = 50;
        mua2 = mua1;
        mua3 = mua1;
        mus2 = mus1;
        mus3 = mus1;
    end

if length(n) == 1
    n1 = n;
    n2 = n;
    n3 = n;
else
    n1 = n(1);
    n2 = n(2);
    n3 = n(3);
end

if isequal(length(time_s),1)
    % DEFAULT VALUES:
    t1 = 0.01e-9; 
    dt = 0.01e-9; % 10 ps time channels
    t2 = 7e-9 + dt; % Up to 7 ns

%     u_end = round((t2-t1)/dt);
%     t(1:u_end) = t1 + (0:u_end-1)*dt;
    time_s = t1:dt:t2-dt;
else
    if time_s(end) > 1
        time_s = time_s / 10^9; % Convert to units seconds
    end
    t1 = time_s(1);
    t2 = time_s(end);
%     dt = time(2)-time(1);
end
u_end = length(time_s);

if ~isempty(varargin) % varargin  allows inputting IRF for convolution
    if length(varargin) == 1
        irf_shifted = varargin{1};
    else
        error('Wrong number of input arguments')
    end
else
    irf_shifted = -1;
end


%% Main body
xh=load('x.dat');
vh=load('w.dat');
Mu = [mua1/n1,mua2/n2,mua3/n3];
del_s = min(Mu);
q1 = 0;
q2 = 30;
qm = (q1+q2)/2;
qr = (q2-q1)/2;
n_max = length(xh);
k_max = 40;
x(1:n_max) = qm - qr*xh(n_max+1-(1:n_max));
q(1:n_max) = qr*vh(n_max+1-(1:n_max));
x((1:n_max)+n_max) = qm + qr*xh(1:n_max);
q((1:n_max)+n_max) = qr*vh(1:n_max);
c0 = 2.99792458*1e11;
del = t2/t1;
alp = 1.09;

A_alp = acosh((pi-2*alp)*del+4*alp-pi)/((4*alp-pi)*sin(alp));
mu = (4*pi*alp-pi^2)/A_alp*k_max/(c0*t2);
h = A_alp/k_max;

[k,n]=meshgrid(0:k_max-1,1:2*n_max);
th_k =(k+1/2)*h;
sk = mu*(1+1i*sinh(th_k+1i*alp));
d_sk = 1i*mu*cosh(th_k+1i*alp);
y = x(n);
w = q(n);

D1 = 1/(3*mus1);
D2 = 1/(3*mus2);
D3 = 1/(3*mus3);
%     D1 = 1/(3*(mus1+mua1)); % Added to test what effect this has
%     D2 = 1/(3*(mus2+mua2));
%     D3 = 1/(3*(mus3+mua3));
z0 = 1/(mua1+mus1);
% zb1 = 2*2.948492562809488*D1; % For n = 1.4
% zb2 = 2*2.948492562809488*D3;

Reff = calc_Reff(n1);
zb1 = (1+Reff)/(1-Reff)*2*D1;
zb2 = (1+Reff)/(1-Reff)*2*D3;


a1 = sqrt(y.^2+(mua1-n1*del_s+n1*sk)/D1);
a2 = sqrt(y.^2+(mua2-n2*del_s+n2*sk)/D2);
a3 = sqrt(y.^2+(mua3-n3*del_s+n3*sk)/D3);

EXP_1 = exp(-a1*(L1+zb1));
EXP_2 = exp(-2*a2*L2);
EXP_3 = exp(-2*a3*(L3+zb2));

T1 = ((n3^2*D2*D3*a2.*a3+(D2*n2*a2).^2).*EXP_3+n3^2*D2*D3*a2.*a3-(D2*n2*a2).^2).*EXP_2+...
    (n3^2*D2*D3*a2.*a3-(D2*n2*a2).^2).*EXP_3+n3^2*D2*D3*a2.*a3+(D2*n2*a2).^2;

T2 = ((n1^2*n3^2*D3*a3+(n1*n2).^2*D2*a2).*EXP_3+n1^2*n3^2*D3*a3-(n1*n2).^2*D2*a2).*EXP_2+...
    (n1^2*n2^2*D2*a2-(n1*n3).^2*D3*a3).*EXP_3-n1^2*n3^2*D3*a3-(n1*n2).^2*D2*a2;

b1 = -exp(-a1*(z0+zb1))./(2*D1*a1);
b2 = -n2^2*exp(-a1*(L1-z0))./(2*D1*a1);
b3 = 1/2*exp(-a1*(L1-z0));

DET    = EXP_1.^2.*(n2^2*T1+D1*a1.*T2) + D1*a1.*T2 - n2^2*T1;
DET_A1 = b1.*(n2^2*T1+D1*a1.*T2).*EXP_1+b3.*T2-b2.*T1;
DET_B1 = EXP_1.*(b2.*T1-b3.*T2)-b1.*(n2^2*T1-D1*a1.*T2);

A1 = DET_A1./DET;
B1 = DET_B1./DET;

R_qs = D1*A1.*a1.*exp(-a1*L1) - D1*B1.*a1.*exp(-a1*zb1);

R_s(1:k_max) = 0;
R_s(:) = 1/(2*pi)*sum(R_qs.*besselj(0,rho*y).*y.*w);

R = zeros(u_end,1);
for u = 1:u_end
    H1 = z0/2*(4*pi*D1*c0/n1)^(-3/2)*time_s(u).^(-5/2).*exp(-mua1*c0/n1*time_s(u));
    H2 = exp(-del_s*c0*time_s(u)+sk(1,:)*c0*time_s(u)).*d_sk(1,:);
    R(u) = c0*h/pi*imag(sum(transpose(R_s(:)).*H2)) + ...
        H1*exp(-(rho^2+z0^2)/(4*D1*c0*time_s(u)/n1));
end


time_ns = time_s * 1e9; % Return time in units of nanoseconds


%% Convolve with IRF
if length(irf_shifted) > 1
    do_plot = 0; % Set this to 1 here to plot DTOF, IRF, and convolved DTOF

    if do_plot == 1
        figure(5); clf
        semilogy(time_ns, R / max(R),'.','Color','red')
        hold on
    end

    R = DTOF_convolve(R, irf_shifted);

    if do_plot == 1
        semilogy(time_ns, R / max(R),'.','Color','blue')
        semilogy(time_ns, irf_shifted / max(irf_shifted),'o','Color','green')
    end
end


end


function R_eff = calc_Reff(n_media) % The angle-averaged probability for reflection at the boundary

    % The used equations are from this Publication:
    % "Boundary conditions for the diffusion equation in radiative transfer"
    % R. C. Haskell, L. O. Svaasand, T. T. Tsay, T. C. Feng, M. McAdams, and B. J. Tromberg,
    % Optical Society of America, A 11, 2727–2741 (1994).
    % DOI:  https://doi.org/10.1364/JOSAA.11.002727


    n_air = 1; % 1 for air, assuming the outside is air

    % Equation 2.3.5
    fun = @(angle) 2 .* sin(angle) .* cos(angle) .* R_fresnel(angle,n_media,n_air);
    R_phi = integral(@(angle)fun(angle), 0, pi/2);

    % Equation 2.3.5
    fun = @(angle) 3 .* sin(angle) .* (cos(angle).^2) .* R_fresnel(angle,n_media,n_air);
    R_j = integral(@(angle)fun(angle), 0, pi/2);

    % Equation 2.3.7
    R_eff = (R_phi + R_j) / (2 - R_phi + R_j);
end


function R_fres = R_fresnel(angle_incident,n_media,n_air) % Equation 2.3.2

    angle_critical = asin( n_air / n_media ); % Using this relation: n_media * sin( angle_critical ) = m_outside

    R_fres = zeros(size(angle_incident,1),1);
    for ind = 1:length(angle_incident)
        % Equation 2.3.2
        if angle_incident(ind) >= angle_critical && angle_incident(ind) <= pi/2
            R_fres(ind) = 1;
    
        else
            % R_s and R_p can be calculated using one of two sets of equations, which are explained on the website:
            % https://en.wikipedia.org/wiki/Fresnel_equations
            % The equations are identical and come from the rules of trigonometry

%             % One way:
%             angle_refracted = asin( n_media/n_air * sin(angle_incident(ind)) ); % Using this relation: n_media sin(angle) = n_outside sin(angle_refracted)
% 
%             % First part of Equation 2.3.2
%             R_s = (n_media*cos(angle_incident(ind)) - n_air*cos(angle_refracted)) / (n_media*cos(angle_incident(ind)) + n_air*cos(angle_refracted));
%             R_s = R_s.^2;
% 
%             % Second part of Equation 2.3.2
%             R_p = (n_media*cos(angle_refracted) - n_air*cos(angle_incident(ind))) / (n_media*cos(angle_refracted) + n_air*cos(angle_incident(ind)));
%             R_p = R_p.^2;
            
            % Second way: This avoids the need to calculate angle_refracted
            val_sqrt = sqrt(1 - (n_media/n_air*sin(angle_incident(ind))).^2);

            % First part of Equation 2.3.2
            R_s = (n_media*cos(angle_incident(ind)) - n_air * val_sqrt ) / (n_media*cos(angle_incident(ind)) + n_air * val_sqrt );
            R_s = abs(R_s);
            R_s = R_s.^2;

            % Second part of Equation 2.3.2
            R_p = (n_media * val_sqrt - n_air*cos(angle_incident(ind))) / (n_media * val_sqrt + n_air*cos(angle_incident(ind)));
            R_p = abs(R_p);
            R_p = R_p.^2;

            R_fres(ind) = R_s / 2 + R_p / 2; % Equation 2.3.2
        end
    end
end