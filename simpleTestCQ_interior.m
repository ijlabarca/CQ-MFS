%% Convolution quadrature:

clc
clear 
close all

% Wave speeds
c1 =1; % interior domain

N = 100;
Np = 20;


M  = 400; % number of times steps
T  = 10;  % final time
dt = T/M; % time increment
tt=(0:dt:T);  
lambda = max([dt^(3/M) eps^(1/2/M)]); % radious complex 

zN = exp(2i*pi/(M+1));
gm = @(z) 0.5*(z.^2-4*z+3);
sl  = gm(zN.^(0:-1:-M)*lambda)/dt;
kl = 1i*sl; % complex frequencies


% We just need to compute half of the frequencies
k_hlf = [kl(1) kl(end:-1:(end-1)/2+2)];


%% Create the geometry
R = 1; Rp = 1.25;
t = linspace(0, 2*pi, N).';

tp = linspace(0, 2*pi, Np).';
% x = R*[cos(t) sin(t)];

% a = 0.6; b = 0.9;
% z = exp(1i*t)- a./(exp(1i*t+b));
% zp = exp(1i*tp)- a./(exp(1i*tp+b));


% a = 0.1; b = 0.9;
% z = R*(exp(1i*t)- a./(exp(1i*t)+b));
% 
% c = 0.4; d = 0.7;
% zp = Rp*(exp(1i*tp)- c./(exp(1i*tp)+d));
% 
% x = [real(z) imag(z)];
% xp =[real(zp) imag(zp)];


% plot(z);hold on
% plot(zp)
% return
% xp = Rp*[cos(tp) 2*sin(tp)];


j = 1;
x = R*[(1+0.3*cos(j*t)).*cos(t) (1+0.3*cos(j*t)).*sin(t)];
xp = Rp*[(1+0.3*cos(j*tp)).*cos(tp) (1+0.3*cos(j*tp)).*sin(tp)];

% verts = [1 1;-1 1;-1 -1;1 -1]; % square
% verts = [0 1;-1 -1;1 -1]; % triangle
% 
% 
% 
% verts = [1 1;-1 0.8;-0.7 -1;1.1 -1;0.9 -0.8]; % square
% 
% verts = [1.5 1.5;0 0.5; -1.5 1.5;-0.5 0;-1.5 -1.5;0 -0.5;1.5 -1.5;0.5 0];
% 
% x = polygon(t, verts); N = (size(verts, 1))*(N-1);
% xp = polygon(tp, 1.25*verts);Np = (size(verts, 1))*(Np-1);


x1 = repmat(x(:, 1),1, Np);
x2 = repmat(x(:, 2),1, Np);

y1 = repmat(xp(:, 1).',N, 1);
y2 = repmat(xp(:, 2).',N, 1);


%% Generate boundary condition

g = zeros(N,M+1); 
dg = zeros(N,M+1);

for n=1:M+1
%     
%     [g(:,n),grad_g] = incident_field(obs.x,tt(n));
%     dg(:,n) = grad_g(:,1).*obs.normals(:,1)+grad_g(:,2).*obs.normals(:,2);
        

    [g(:,n),~] = incident_field(x,tt(n));
%     g(:, n) = cos(3*pi*(x(:,1) + x(:, 2)) - c1*tt(n)).*(1-POU(tt(n), 1, T/2))*(tt(n)>0);
%     dg(:,n) = grad_g(:,1).*obs.normals(:,1)+grad_g(:,2).*obs.normals(:,2);

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

    k1 = k_hlf(n)/c1;
    
    A = 1i/4*besselh(0, k1*sqrt((x1-y1).^2 + (x2-y2).^2));

    phip_hlf(:,n) = A\g_hlf(:,n);
    disp(n);
    
end

%% Evaluate the field

Xlim = [-2 2];
Ylim = [-2 2];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*20));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*20));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));


up_hlf = zeros(size(pts,1),M/2+1);


for n=1:M/2+1
    
    k1 = k_hlf(n)/c1;  
    for m = 1:Np
    up_hlf(:, n) = up_hlf(:, n)+ phip_hlf(m, n)*1i/4*besselh(0, k1*sqrt((pts(:,1)-xp(m, 1)).^2 + (pts(:,2)-xp(m, 2)).^2));
    end
    disp(n);

end

up = [up_hlf(:,1) conj(up_hlf(:,2:end)) up_hlf(:,end:-1:2)];


% Inverting Z-trasnform
Lam = repmat(lambda.^(0:M),size(pts,1),1);

%%

% uinc = zeros(size(pts,1),M+1); 
% for n = 1:M+1
% [uinc(:, n), ~] = incident_field(pts,tt(n));
% end
u = Lam.^(-1).*ifft(up,[],2);
% u(in, :) = 1;
u(~in, :) = 0;

%% Plot the solution
close all  

X = reshape(pts(:,1),Ny,Nx);

Y = reshape(pts(:,2),Ny,Nx);


for n=1:M+1
    
    
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx),'EdgeColor','none');hold on
    plot3([x(:,1); x(1,1)],[x(:,2); x(1,2)],1*ones(N+1,1),'LineWidth',2,'color','k');
    hold off
    view(2)
    shading interp;
    title(n)
%     colormap(brewermap([],'*RdBu'))
    colormap hot;
    caxis([-1 1])
    colorbar
    axis equal
    axis tight
    
    drawnow
    
    pause(dt)

end



