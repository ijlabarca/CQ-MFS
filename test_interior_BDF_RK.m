%% Method of Fundamental Solution + Convolution Quadrature
%  Time domain test for the interior problem
%  with exact solutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ignacio Labarca 
% Seminar for Applied Mathematics, ETH Zurich
% email:  ignacio.labarca@sam.math.ethz.ch
% date:   July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% T      : Final time
% M      : Number of time steps
% dt     : time step
% bdf    : BDF1 (1) or BDF2 (2) formulas for multistep CQ
% RK     : RadauIIA of two stages (1) or three stages (2) for multistage CQ
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
tic 
warning('off');


%%%%%%%%%%%%%%%%%%%%%%
% Problems:
%%%%%%%%%%%%%%%%%%%%%%
% Circle           : 0
% Rounded Triangle : 1
% Inverted Ellipse : 2
%%%%%%%%%%%%%%%%%%%%%%
problem = 0;


% Wave speed
wavespeed =1; % exterior domain

Ms = [50 100 200 400 800 1600];

R = 1;     
Rp = 1.2;  

a1 = 0.3;
a2 = 0.25;

if problem == 0
    Z = @(z) z; % Circle
elseif problem == 1
    Z = @(z) z+a1./(z.^2); % Rounded triangle;
else
    Z = @(z) z./(1+a2.*z.^2); % Inverted ellipse;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF NUMERICAL EXPERIMENTS

T  = 10;  % Time interval [0, T]

for bdf = 1:2
count = 1;
for M = Ms
    
    
dt = T/M; % time increment
tt=(0:dt:T);  
lambda = eps^(1/2/M);

zN = exp(2i*pi/(M+1));

if bdf > 1
gm = @(z) 0.5*(z.^2-4*z+3);      %BDF2 
else
gm = @(z) 1-z;                   %BDF1
end

sl  = gm(zN.^(0:-1:-M)*lambda)/dt;
kl = 1i*sl; % complex wavenumbers

% We compute half of the wavenumbers
k_hlf = [kl(1) kl(end:-1:(end-1)/2+2)];

Np = 500; 
N = 2*Np;
%% Create the geometry

expN = exp(2i*pi/N);
z1 = Z(R*expN.^(0:N-1)).';
x = [real(z1) imag(z1)];

expNp = exp(2i*pi/Np);
z2 = Z(Rp*expNp.^(0:Np-1)).';
xp = [real(z2) imag(z2)];

%% Plot Geometry
%%%%%%%%%%%%%%%%%%%

scatter(x(:, 1), x(:, 2), 'filled', 'ko'); hold on;
scatter(xp(:, 1), xp(:, 2), 'filled','rs'); axis off;hold off;
pause(5);
close all

%%%%%%%%%%%%%%%%%%%
x1 = repmat(x(:, 1),1, Np);
x2 = repmat(x(:, 2),1, Np);

y1 = repmat(xp(:, 1).',N, 1);
y2 = repmat(xp(:, 2).',N, 1);

% Generate boundary condition
g = zeros(N,M+1); 
for n=1:M+1
    [g(:,n),~] = incident_field(x,tt(n));
end

g = -g;
% Computing the Z-transform of boundary data
Lam = repmat(lambda.^(0:M),N,1);

gp = fft(Lam.*g,[],2);

% again, we just need hald of the frequencies
g_hlf=[gp(:,1) gp(:,end:-1:(end-1)/2+2)];

%% Frequency domain problems
alphap_hlf = zeros(Np,M/2+1);

hbar = parfor_progressbar(M,['Solving linear systems... M = ',num2str(M/2)]);
for n=1:M/2+1

    k = k_hlf(n)/wavespeed;
    A = 1i/4*besselh(0, k*sqrt((x1-y1).^2 + (x2-y2).^2));
    
    if min(min(abs(A))) < 1e-20
        A(abs(A) < 1e-20) = 0;
        A = sparse(A);
    end
    
    alphap_hlf(:,n) = A\g_hlf(:,n);
    hbar.iterate(1);
    
end
close(hbar);

clc
%% Evaluate the field

% We evaluate the solution 'up' in the set of points 'pts'.
% then, we recover time domain solution 'u'.



if problem == 0 || problem == 1
    pts = 0.25*[2 2;2 -2;-2 -2;-2 2];
elseif problem == 1
    pts = 0.25*[-2 -2;-2 2;0 0;2 0];
end

Ntot = size(pts,1);

in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));


up_hlf = zeros(Ntot,M/2+1);


for n=1:M/2+1
    k = k_hlf(n)/wavespeed;  
    for m = 1:Np
    up_hlf(:, n) = up_hlf(:, n)+ alphap_hlf(m, n)*1i/4*besselh(0, k*sqrt((pts(:,1)-xp(m, 1)).^2 + (pts(:,2)-xp(m, 2)).^2));
    end
end

up = [up_hlf(:,1) conj(up_hlf(:,2:end)) up_hlf(:,end:-1:2)];

% Inverse Z-transform
Lam = repmat(lambda.^(0:M),Ntot,1);

u = Lam.^(-1).*ifft(up,[],2);


% Reference solution
uinc = zeros(size(u));
for n=1:M+1 
[uinc(:,n),~] = incident_field(pts,wavespeed*(tt(n)));
end

u(in, :) = u(in, :) + uinc(in, :);
u(~in,:)= 0;

% Compute errors
sum =0;
norm_u = 0;

for m = 1:M+1
   sum = sum + norm(u(:, m))^2;
   norm_u = norm_u + norm(uinc(:, m))^2;
end
error(count, bdf)= sqrt(sum)/sqrt(norm_u)
if count > 1
    order(count, bdf) = -log(error(count,bdf)/error(count-1,bdf))/log(Ms(count)/Ms(count-1));
end

count = count+1;

pause(2)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUNGE - KUTTA NUMERICAL EXPERIMENTS

% RK = 1 : RadauIIA - two stages
% RK = 2 : RadauIIA - three stages

for RK = 1:2
count= 1;

if RK ==1
loop = Ms; 
else
  
loop = [50 100 200 400 800];
end

for M = loop
dt = T/M; % time increment
tt=(0:dt:T);
lambda = eps^(1/2/M);

zN = exp(2i*pi/M);

if RK == 1
A_RK=[5/12 -1/12; 3/4 1/4]; % Radau IIa 2nd order     

else
A_RK=[11/45 37/225 -2/225; ...                     % Radau IIa 5th order
     37/225 11/45 -2/225; ...
     4/9 4/9 1/9] + ...
    [-7*sqrt(6)/360 -169*sqrt(6)/1800 sqrt(6)/75; ...
      169*sqrt(6)/1800 7*sqrt(6)/360 -sqrt(6)/75; ...
     -sqrt(6)/36 sqrt(6)/36 0];

end

b=A_RK(end,:);
stages=size(A_RK,1);
invA = inv(A_RK);
c=A_RK*ones(stages,1);
B = (invA*ones(size(A_RK,1),1))*[zeros(1,stages-1),1];

Np = 500;
N = 2*Np;
%% Create the geometry

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
idx_N=@(s) (s-1)*N+1:s*N;
idx_Np=@(s) (s-1)*Np+1:s*Np;


%% Right Hand side

F =zeros(N,stages,M);
for st=1:stages
    
    for n=1:M
    [g1(:,n),~] = incident_field(x,wavespeed*(tt(n)+dt*c(st)));
    end
    F(:,st, :) = g1; clear g1
end
F = reshape(F, [stages*N, M]);

%% Computing the Z-transform of boundary data
Lam = repmat(lambda.^(0:M-1),stages*N,1);
F = fft(Lam.*F,[],2);

%% Solution of the problem in the frequency domain
alphap_hlf = zeros(Np*stages,M);

hbar = parfor_progressbar(floor(M/2),['Solving linear systems... M = ',num2str(M/2)]);
tic
for n=0:floor(M/2)
    [P,Lambda]=eig(invA-lambda*zN^(-(n))*B);
    Lambda=diag(Lambda)/dt;
    gl=kron(inv(P),speye(N))*F(:,n+1);

    ul=zeros(stages*Np,1);

    for s=1:stages
        k = 1i*Lambda(s)/wavespeed;
        A = 1i/4*besselh(0, k*sqrt((x1-y1).^2 + (x2-y2).^2));
        if min(min(abs(A))) < 1e-20
            A(abs(A) < 1e-20) = 0;
            A = sparse(A);
        end
        ul(idx_Np(s))=A\gl(idx_N(s));
    end

    alphap_hlf(:,n+1) = ul;

hbar.iterate(1);

end
toc
close(hbar);

alphap_hlf(:,M+1-(1:floor((M-1)/2))) = conj(alphap_hlf(:,2:floor((M-1)/2)+1));
%% Evaluate the field



if problem == 0 || problem == 1
    pts = 0.25*[2 2;2 -2;-2 -2;-2 2];
elseif problem == 1
    pts = 0.25*[-2 -2;-2 2;0 0;2 0];
end

Ntot = size(pts, 1);

in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));

idx=@(s) (s-1)*Ntot+1:(s*Ntot); 
idy=@(s) (s-1)*Np+1:(s*Np);

u_hlf = zeros(stages*Ntot,M);
tic
for n=1:M
    [P,Lambda]=eig(invA-lambda*zN^(-(n-1))*B);
    Lambda=diag(Lambda)/dt;

    gl=alphap_hlf(:,n);
    ul=zeros(stages*Ntot,1);

    for s=1:stages
        
    u_s = zeros(Ntot, 1);
    k = 1i*Lambda(s)/wavespeed;
    alpha = gl(idy(s));
    for m = 1:Np
        u_s = u_s+ alpha(m)*1i/4*besselh(0, k*sqrt((pts(:,1)-xp(m, 1)).^2 + (pts(:,2)-xp(m, 2)).^2));
    end
    ul(idx(s))=u_s;
    end
    u_hlf(:,n)=kron(P,speye(Ntot))*ul;

end

% Inverting Z-transform
Lam = repmat(lambda.^(0:M-1),stages*Ntot,1);


u = Lam.^(-1).*ifft(u_hlf,[],2);
u = u((stages-1)*Ntot+1:end,:);
u = [zeros(size(u, 1), 1) u];

uinc = zeros(size(u));
for n=2:M+1
[uinc(:,n),~] = incident_field(pts,wavespeed*(tt(n)));
  
end

u(in, :) = u(in, :) - uinc(in, :);
u(~in,:)= 0;
%% Compute errors
% L2 norm in [0, T], pointwise error in 'pts'

sum =0;
norm_u = 0;

for m = 1:M+1
   sum = sum + norm(u(:, m))^2;
   norm_u = norm_u + norm(uinc(:, m))^2;
end
error(count, RK+bdf)= sqrt(sum)/sqrt(norm_u)
if count > 1
    order(count, RK+bdf) = -log(error(count,RK+bdf)/error(count-1,RK+bdf))/log(Ms(count)/Ms(count-1));
end
count = count+1;
pause(3)
end
toc
end

%% Convergence Plot

loglog(Ms, error(:,1), '--gx', ...
      Ms, error(:,2), '--md', ...
      Ms, error(:,3), '--ro', 'Linewidth',2, 'MarkerSize', 14); hold on;
loglog(Ms, error(:,4), '--bs', 'Linewidth',2, 'MarkerSize', 14);


loglog(Ms, 2.2e+2*Ms.^(-1), '-k', ...
      Ms, 0.2e+4*Ms.^(-2), '-k', ...
      Ms, 2.2e+4*Ms.^(-3), '-k', ...
      Ms, 0.2e+7*Ms.^(-5), '-k', 'Linewidth', 2);
text(800, 0.5, '1','FontSize', 18);

text(700, 0.0020, '2','FontSize', 18);

text(600, 0.00022, '3','FontSize', 18);

text(300, 0.000002, '5','FontSize', 18);
set(gca,'FontSize', 14);
legend({'BDF1','BDF2','RadauIIA - two stages', 'RadauIIA - three stages'},'Location','southwest', 'FontSize', 16);
xlabel('Number of time steps','Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative error in $L^2-$norm', 'Interpreter', 'latex', 'FontSize', 18)
xticks(Ms)
ylim([0.2e-9 10]);
xlim([0 1600]);
title({'Convergence of interior', 'problem'},'Interpreter', 'latex', 'FontSize', 24)