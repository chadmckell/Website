%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Lofelt - Membrane Filter Animation
% 
% Author: Chad McKell
% Date: 28 June 2017
%
% Description: Animation of a physically-modeled membrane filter for an  
% audio input signal. The modal synthesis algorithm includes forcing, 
% frequency-dependent loss, and zero stiffness. The forcing is provided by 
% the input audio signal F. The animation is saved as 'membrane.avi'.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear; close all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define forcing term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Read audio file
[F,fs] = audioread('shot.wav'); % gun, kick, shot, square, tom
% F = zeros(100000,1);
% F(100) = 1;
% fs = 44100;

% If file contains stereo signal, convert to mono
if size(F,2) > 1
    nz = ~sum(F == 0,2); % Find rows with non-zero values
    F = sum(F,2); % Add columns of F 
    F(nz) = F(nz)/2; % Average non-zero rows
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Set exposable variables 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Membrane parameters 
R = 1; % aspect ratio (1 to 10)
c = 67; % fundamental frequency (15 to 1000)

% Algorithm parameters
N = 18; % number of modes (1 to 10)
dR = 0.4; % decay rate (0.5 to 5)
dS = 2; % decay scale (0.5 to 10)
lxi = 0.2; % excitation location ratio along x (0.01 to 1)
lyi = 0.6; % excitation location ratio along y (0.01 to 1)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute derived variables 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
p = (1:N)'; % mode grid along x
q = (1:N)'; % mode grid along y
Lx = 1; % x-coordinate length [meter] 
Ly = Lx/R; % y-coordinate length [meter]
env = exp(-dR*p*q'); % nonlinear decay envelope
T60 = env*dS; % frequency-dependent reverberation times [s]
Nf = length(F); % length of input signal [samples]
k = 1/fs; % sample duration [sec]
sig = 6*log(10)./T60; % damping coefficients
xi = lxi*Lx; % excitation location along x [meter]
yi = lyi*Ly; % excitation location along y [meter]

% Compute modal frequencies of surface
w = zeros(N,N);
for n=1:N
    for m=1:N
        w(n,m) = sqrt(n.^2/Lx^2 + m.^2/Ly^2)*pi*c; 
    end
end
f0 = w(1,1)/(2*pi);

xgrid = (1:90)/90; % vector of grid points for plotting
ygrid = (1:90)/90; % vector of grid points for plotting

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define coefficients in SHO difference scheme
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
frac = 1./(1+sig*k);
a1 = frac.*(2-w.^2*k^2);
a2 = frac.*(sig*k - 1);
a3 = k*sin(p*pi*xi/Lx)*sin(q'*pi*(yi)/Ly);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize harmonic oscillator vectors 
phi = zeros(N,N); % current time step
phi1 = zeros(N,N); % 1 time step back
phi2 = zeros(N,N); % 2 time steps back

% Set variable to saving video
v = VideoWriter('membrane.avi');

Ux = sin(p.*pi*xgrid/Lx)';
Uy = sin(q.*pi*ygrid/Ly);

tic
for n = 1:Nf
    
    % Compute time-dependent harmonic oscillator component for each mode
    phi = a1.*phi1 + a2.*phi2 + a3*F(n);

    % Calcuate output signal then plot the signal
    outP = Ux*phi*Uy;

    if ~mod(n,5)
        surf(xgrid, ygrid, outP)
        zlim([-8*10^-1 8*10^-1]);
        time = n*k;
%         title(['t = '  num2str(time) ' sec'])
        axis off
        pbaspect([1 1 0.81])
        drawnow
        open(v);
        writeVideo(v, getframe(gca));
    end
    
    % Set harmonic oscillator vectors equal to next grid line in time
    phi2 = phi1; 
    phi1 = phi;
end
toc

close(v)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Bilbao, S. "Numerical Sound Synthesis: Finite Difference Schemes and
%    Simulation in Musical Acoustics". John Wiley & Sons (2009).