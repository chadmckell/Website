%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Lofelt - String Filter
% 
% Author: Chad McKell
% Date: 22 May 17
%
% Description: A physically-modeled string filter for an audio input signal.
% The filter is a modal-synthesized string with forcing, frequency-
% dependent loss, and zero stiffness. The forcing is provided by the input 
% signal F. 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear; close all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define forcing term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Read input signal (single impulse)
[F,fs] = audioread('shot.wav'); % gun, kick, shot, square, tom
% F = zeros(100000,1);
% F(10) = 1;
% fs = 44100;

% If file contains stereo signal, convert to mono
if size(F,2) > 1
    nz = ~sum(F == 0,2); % Find rows with non-zero values
    F = sum(F,2); % Add columns of F 
    F(nz) = F(nz)/2; % Average non-zero rows
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Set global variables 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% String parameters 
L = 0.4; % length of string [meter] 
Ts = 11; % tension [N]
r = 0.00034; % string radius [m]
rho = 6504; % density [kg/m^3]

% Algorithm parameters
N = 9; % number of modes
dR = -0.01; % decay rate (-1.5 to 0.1)
dS = 10; % decay scale (0.5 to 10)
li = 0.134; % xi ratio (0.01 to 1) 
lo = 0.5; % xo ratio (0.01 to 1) 
Np = 100;  % length of spatial grid points for plotting
u0 = 0; % initial displacement of string [meter]
v0 = 4; % inital velocity of string [meter/sec]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Triangular excitation function 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% Initialize excitation vector
rc = zeros(N,1);

% triangular excitation
% A = 1;
% for q = 1:N
%     x = L*(q-1)/N;
%     if x <= L/3
%         rc(q) = (3/10)*x/L;
%     else
%         rc(q) = (3*(1-(x/L))/20);
%     end
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute derived variables 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
m = (1:N)'; % mode grid
y = exp(dR*m); % nonlinear decay envelope
T60 = y*dS; % frequency-dependent reverberation times [s]
Nf = length(F); % length of input signal [samples]
k = 1/fs; % sample duration [sec]
As = pi*r^2; % cross-sectional area of string
c = 9;%sqrt(Ts/(rho*As)); % wave speed [meter/sec]
f0 = c/(2*L); % fundamental frequency of vibration [Hz]
sig = 6*log(10)./T60; % damping coefficients
xgrid = (0:Np)/Np; % vector of grid points for plotting
xi = li*L; % excitation location [meter]
xo = lo*L; % audio read out location [meter]
wp = m*2*pi*f0; % angular frequencies of modes
fp = m*f0; % frequencies of modes

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define coefficients in SHO difference scheme
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
frac = 1./(1+sig*k);
a1 = frac.*(2-wp.^2*k^2);
a2 = frac.*(sig*k - 1);
a3 = k*sin(m*pi*xi/L);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check stability (poles must lie inside the unit circle) 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
poles = a1/2 + sqrt((a1.^2+4*a2)/4);
stab = sqrt(real(poles).^2 + imag(poles).^2);
for n=1:length(stab)
    if stab(n) > 1
        error('System is unstable. Check poles.')
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute DC gain
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
gain = a3./(1-a1-a2);
tot = sum(gain);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize output signals
out = zeros(Nf,1); % audio signal

% Initialize harmonic oscillator vectors 
phi = zeros(N,1); % current time step
phi1 = zeros(N,1); % 1 time step back
phi2 = zeros(N,1); % 2 time steps back

% Set spatial vector
Up = sin(m.*pi*xgrid); % plot vector

% Compute coefficients of SHO analytical form (see Ref 1, p. 153)
rcfs = -imag(fft([rc; zeros(N,1)])); 
rcfs = 2*rcfs(2:N+1)/N; 

% Set initial conditions
phi1 = (u0*cos(wp*k) + v0*sin(wp*k)./wp).*rcfs;
phi2 = u0*rcfs;

% Set variable to saving video
v = VideoWriter('string.avi');

total1 = [];
total2 = [];
total3 = [];
total4 = [];
total5 = [];
total6 = [];

tic
for n = 1:Nf
    
    % Compute time-dependent harmonic oscillator component for each mode
    phi = a1.*phi1 + a2.*phi2 + a3.*F(n);
    
    % Calculate output signal at xo 
    out = Up.*phi;
    
    if ~mod(n,10)
        total = sum(out);
        time = n*k;
        plot(xgrid, total, 'k', 'LineWidth',1)
        ylim([-3*10^0 3*10^0])
        axis off
        pbaspect([1 1 0.81])
%         title(['t = '  num2str(time) ' sec'])
    %     xlabel('Length (normalized)')
    %     ylabel('Amplitude')
        drawnow
        open(v);
        writeVideo(v, getframe(gca));
    end

%     if n == 1
%         new=0.2526*sum(out(1:40,:));
%         total1 = [total1; new];
%     end
%     if n == 1000
%         new=0.2526*sum(out(1:40,:));
%         total2 = [total2; new];
%     end
%     if n == 2000
%         new=0.2526*sum(out(1:40,:));
%         total3 = [total3; new];
%     end
%     if n == 3000
%         new=0.2526*sum(out(1:40,:));
%         total4 = [total4; new];
%     end
%     if n == 4000
%         new=0.2526*sum(out(1:40,:));
%         total5 = [total5; new];
%     end
%     if n == 5000
%         new=0.2526*sum(out(1:40,:));
%         total6 = [total6; new];
%     end
    
    % Set harmonic oscillator vectors equal to next grid line in time
    phi2 = phi1; 
    phi1 = phi;
end
toc


% close(v);



% Basic_Wave_Equation_Script.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SR2 = 44100;   % sampling rate
c2 = c;        % speed of wave
L2 = 1;        % length of string
lambda2 = 1;     % courant number
v02 = 0;
% N=390; h = L/N; lambda = c*k/h;

%%%% flags
excite2 = 0; % 1: hann; 0: triangular
bc2 = 0; % if 0, Dirichlet; if 1, Neumann.
%%%%

xc2 = .1; % striking point along x?
xout2 = .5; % point to read out along string
A2 = .01;
wid2 = .02; % hann width

Tf2 = 1; % simulation time
k2 = 1/SR2; % time step
Nf2 = floor(Tf2*SR2); % number of samples

% Need to have an integer multiple N grid points
h2 = c2*k2/lambda2; % grid spacing
N2 = floor(L2/h2); % Number of grid points
h2 = L2/N2; % redefine h
lambda2 = c2*k2/h2; % redefine lambda


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II Generate Hann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rc2 = zeros(N+1,1);

if excite2
    for q = 1:N2+1
      x = (q-1)*h2;
      dist = abs(x-xc2);
      if dist <= wid2
        rc2(q) = .5*A2*(1+cos(pi*dist/wid2));
      end
    end
else
   for q = 1:N2+1
       x = (q-1)*h2;
       if x <= L2/3
           rc2(q) = (3/10)*x;
       else
           rc2(q) = 3*(1-x)/20;
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III Initialise u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u2 = rc2;
u1 = rc2; % initial velocity is zero in this case
u = zeros(N2+1,1);

out2 = u;

tic
for n = 1:1
    
  if bc2
    u(1) = (2-lambda2^2)*u1(1) + lambda2^2*u1(2) - u2(1);
    u(N2+1) = (2-lambda2^2)*u1(N2+1) + lambda2^2*u1(N2-1) - u2(N2+1);
  end

  % this is the difference equation pretty much straight from the book.
  u(2:N2) = 2*(1-lambda2^2)*u1(2:N2) + ((lambda2^2)*(u1(3:N2+1) + u1(1:N2-1))) - u2(2:N2);
  out2(n) = u(floor(xout2*N2/L2));
  
  u2 = u1;
  u1 = u;
end
toc

total3m = 0.2526*sum(out(1:3,:));
total6m = 0.2526*sum(out(1:6,:));
total40m = 0.2526*sum(out(1:40,:));
out(1,:) = 0.2526*out(1,:);
out(2,:) = 0.2526*out(2,:);
out(3,:) = 0.2526*out(3,:);
out(4,:) = 0.2526*out(4,:);
out(5,:) = 0.2526*out(5,:);
out(6,:) = 0.2526*out(6,:);

% out(1,:) = out(1,:)/max(out(1,:));
% out(2,:) = out(2,:)/max(out(2,:));
% out(3,:) = out(3,:)/max(out(3,:));
% out(4,:) = out(4,:)/max(out(4,:));
% out(5,:) = out(5,:)/max(out(5,:));
% out(6,:) = out(6,:)/max(out(6,:));

% fig = figure(1);
% 
% subplot_tight(3,2,1,[0.1])
% plot(xgrid, out(1,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-1 0 1])
% ylim([-1 1])
% title('p = 1');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,2,[0.1])
% plot(xgrid, out(4,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% title('p = 4');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,3,[0.1])
% plot(xgrid, out(2,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% title('p = 2');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,4,[0.1])
% plot(xgrid, out(5,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% title('p = 5');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,5,[0.1])
% plot(xgrid, out(3,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% title('p = 3');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,6,[0.1])
% plot(xgrid, out(6,:), 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% title('p = 6');
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% tightfig(fig);
% saveas(fig,'z','epsc')
% 
% 
% 
% 
% grid = [0:N]/N;
% grid = grid(1:end-1);
% 
% fig2 = figure(2);
% % plot([0:N2]'*h2,u,'k', 'LineWidth',1)
% % hold on
% plot(xgrid, total3m, ':b','LineWidth',2)
% hold on
% plot(xgrid, total6m, '--r', 'LineWidth',2)
% hold on
% plot(xgrid, total40m, 'k', 'LineWidth',2) %'Color', [0.6 0.6 0.6]
% % hold on
% % plot(grid, rc,'k', 'LineWidth',2)
% % xticks([0 0.5 1])
% % yticks([-0.1 0 0.1])
% % ylim([-0.1 0.1])
% % title('Sum (50 modes)');
% set(gca,'fontsize',16)
% % pbaspect([2 1 1])
% 
% xticks([0 0.5 1])
% yticks([0 0.1])
% ylim([0 0.11])
% title('Sums');
% set(gca,'fontsize',16)
% legend('3 modes', '6 modes', '40 modes');
% pbaspect([2 1 1])
% 
% tightfig(fig2);
% saveas(fig2,'sum','epsc')
% 
% exact = figure(3)
% plot(grid, rc,'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([0 0.1])
% ylim([0 0.11])
% xlabel('x')
% ylabel('g(x)')
% set(gca,'fontsize',16)
% pbaspect([2 1 1])
% tightfig(exact);
% saveas(exact,'exact','epsc')
% 
% 
% 
% 
% time = figure(4);
% subplot_tight(3,2,1,[0.1])
% plot(xgrid, total1, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% title(['t = 0.00 sec' ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,3,[0.1])
% plot(xgrid, total2, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% nk = 1000*k;
% title(['t = 0.02 sec']);%  num2str(nk) ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,5,[0.1])
% plot(xgrid, total3, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% nk = 2000*k;
% title(['t = 0.05 sec']);%'  num2str(nk) ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,2,[0.1])
% plot(xgrid, total4, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% nk = 3000*k;
% title(['t = 0.07 sec']);%'  num2str(nk) ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,4,[0.1])
% plot(xgrid, total5, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% nk = 4000*k;
% title(['t = 0.09 sec']);%'  num2str(nk) ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% subplot_tight(3,2,6,[0.1])
% plot(xgrid, total6, 'k', 'LineWidth',2)
% xticks([0 0.5 1])
% yticks([-0.1 0 0.1])
% ylim([-0.11 0.11])
% nk = 5000*k;
% title(['t = 0.11 sec']);%'  num2str(nk) ]);
% set(gca,'fontsize',16)
% pbaspect([3 1 1])
% 
% tightfig(time);
% saveas(time,'time','epsc')



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Bilbao, S. "Numerical Sound Synthesis: Finite Difference Schemes and
%    Simulation in Musical Acoustics". John Wiley & Sons (2009).
