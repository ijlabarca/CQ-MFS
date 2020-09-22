%% Convolution quadrature + Method of Fundamental Solutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ignacio Labarca 
% Seminar for Applied Mathematics, ETH Zurich
% email:  ignacio.labarca@sam.math.ethz.ch
% date:   July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c      : wavespeed
% T      : Final time
% M      : Number of time steps
% dt     : time step
% R      : Radius for collocation points in Gamma 
% Rp     : Radius for charge points in Sigma
% Z      : Parametrization of boundary Gamma
% lambda : Radius of circle from contour integral
% gm     : Rational polynomial for multistep CQ
% x      : Collocation points
% xp     : Charge points



clc
clear 
close all




%%%%%%%%%%%%%%%%%%%%%%
% Problems:
%%%%%%%%%%%%%%%%%%%%%%
% Circle           : 0
% Rounded Triangle : 1
% Inverted Ellipse : 2
%%%%%%%%%%%%%%%%%%%%%%
problem = 0;




% Wave speeds
wavespeed =1; % exterior domain

N = 600;
Np = 300;

M  = 200; % number of times steps
Mmax = M;
T  = 10;  % final time
dt = T/M; % time increment
tt=(0:dt:T);  
lambda = eps^(1/2/M); % radious complex 

zN = exp(2i*pi/(M+1));
gm = @(z) 0.5*(z.^2-4*z+3);
sl  = gm(zN.^(0:-1:-M)*lambda)/dt;
kl = 1i*sl; % complex frequencies


% We just need to compute half of the frequencies
k_hlf = [kl(1) kl(end:-1:(end-1)/2+2)];


%% Create the geometry
R = 1; Rp = 0.8;
t = linspace(0, 2*pi, N).';


tp = linspace(0, 2*pi, Np).';



a1 = 0.3;
a2 = 0.25;

if problem == 0
    Z = @(z) z; % Circle
elseif problem == 1
    Z = @(z) z+a1./(z.^2); % Rounded triangle;
else
    Z = @(z) z./(1+a2.*z.^2); % Inverted ellipse;
end

expN = exp(2i*pi/N);
z1 = Z(R*expN.^(0:N-1)).';
x = [real(z1) imag(z1)];


expNp = exp(2i*pi/Np);
z2 = Z(Rp*expNp.^(0:Np-1)).';
xp = [real(z2) imag(z2)];


x1 = repmat(x(:, 1),1, Np);
x2 = repmat(x(:, 2),1, Np);

y1 = repmat(xp(:, 1).',N, 1);
y2 = repmat(xp(:, 2).',N, 1);


%% Generate boundary condition

g = zeros(N,M+1); 

for n=1:M+1
    [g(:,n),~] = incident_field(x,tt(n));
end

g = -g;
%% Computing the Z-transform of boundary data
Lam = repmat(lambda.^(0:M),N,1);

gp = fft(Lam.*g,[],2);

% again, we just need hald of the frequencies
g_hlf=[gp(:,1) gp(:,end:-1:(end-1)/2+2)];


%% Solution of the problem in the frequency domains 
phip_hlf = zeros(Np,M/2+1);

for n=1:M/2+1

    k1 = k_hlf(n)/wavespeed;
    
    A = 1i/4*besselh(0, k1*sqrt((x1-y1).^2 + (x2-y2).^2));

    phip_hlf(:,n) = A\g_hlf(:,n);
    disp(n);
    
end

%% Evaluate the field

Xlim = [-3 3];
Ylim = [-3 3];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*10));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*10));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));


up_hlf = zeros(size(pts,1),M/2+1);


for n=1:M/2+1
    
    k1 = k_hlf(n)/wavespeed;  
    for m = 1:Np
    up_hlf(:, n) = up_hlf(:, n)+ phip_hlf(m, n)*1i/4*besselh(0, k1*sqrt((pts(:,1)-xp(m, 1)).^2 + (pts(:,2)-xp(m, 2)).^2));
    end
    disp(n);

end

up = [up_hlf(:,1) conj(up_hlf(:,2:end)) up_hlf(:,end:-1:2)];


% Inverting Z-trasnform
Lam = repmat(lambda.^(0:M),size(pts,1),1);

%%

uinc = zeros(size(pts,1),M+1); 
for n = 1:M+1
[uinc(:, n), ~] = incident_field(pts,tt(n));
end

%%
u = Lam.^(-1).*ifft(up,[],2);
u(~in, :) =u(~in, :)+uinc(~in, :); 

u(in, :) = 0;

%% Plot the solution
close all  

X = reshape(pts(:,1),Ny,Nx);

Y = reshape(pts(:,2),Ny,Nx);





for n=1:M+1
    
    
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx),'EdgeColor','none');hold on
    plot3([x(:,1); x(1,1)],[x(:,2); x(1,2)],1*ones(N+1,1),'LineWidth',2,'color','k');
    plot3([xp(:,1); xp(1,1)],[xp(:,2); xp(1,2)],1*ones(Np+1,1),'LineWidth',2,'color','r');
    
    hold off
    view(2)
    shading interp;
    colormap jet;
    caxis([-1 1])
%     colorbar
    axis equal
    axis tight
    
    set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
    set(gca,'Yticklabel',[]) %to just get rid of the numbers but leave the ticks.
    
    drawnow
    
    pause(dt)

end



