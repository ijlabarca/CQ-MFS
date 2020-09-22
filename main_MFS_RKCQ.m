%% Runge-Kutta Convolution quadrature + Method of Fundamental Solutions


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
% x      : Collocation points
% xp     : Charge points



clc
clear 
close all

% Wave speeds
c(1) = 1; % exterior domain


N = 1000;
Np = 500;

M  = 400; % number of times steps
T  = 10;  % final time
dt = T/M; % time increment
tt=(0:dt:T);  
lambda = eps^(1/2/M); % radious complex 
tol = 1e-20;
zN = exp(2i*pi/M);

% Ark=[5/12 -1/12; 3/4 1/4]; % Radau IIa 2nd order


Ark=[11/45 37/225 -2/225; ...                     % Radau IIa 5th order
     37/225 11/45 -2/225; ...
     4/9 4/9 1/9] + ...
    [-7*sqrt(6)/360 -169*sqrt(6)/1800 sqrt(6)/75; ...
      169*sqrt(6)/1800 7*sqrt(6)/360 -sqrt(6)/75; ...
     -sqrt(6)/36 sqrt(6)/36 0]; 
b=Ark(end,:);
S=size(Ark,1);
Am1 = inv(Ark);
crk=Ark*ones(S,1);
B = (Am1*ones(size(Ark,1),1))*[zeros(1,S-1),1];

%% Create the geometry
R = 1; Rp = 0.8;


a1 = 0.3;
a2 = 0.25;
% Z = @(z) z;
% Z = @(z) z+a1./(z.^2); % Rounded triangle;
Z = @(z) z./(1+a2.*z.^2); % Inverted ellipse;

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



idx_N=@(s) (s-1)*N+1:s*N;
idx_Np=@(s) (s-1)*Np+1:s*Np;


%% Right Hand side
F =zeros(N,S,M);

for st=1:S
%     
    for n=1:M
    
    [g1(:,n),~] = incident_field(x,c(1)*(tt(n)+dt*crk(st)));
      
    end

    F(:,st, :) = g1;   
    clear g1
end

F = reshape(F, [S*N, M]);




%% Computing the Z-transform of boundary data
Lam = repmat(lambda.^(0:M-1),S*N,1);
F = fft(Lam.*F,[],2);



%% Solution of the problem in the frequency domains 
phip_hlf = zeros(Np*S,M);


hbar = parfor_progressbar(floor(M/2),['Solving linear systems... M = ',num2str(M/2)]);
tic
for n=0:floor(M/2)
    

    [P,Lambda]=eig(Am1-lambda*zN^(-(n))*B);   
    Lambda=diag(Lambda)/dt;
    gl=kron(inv(P),speye(N))*F(:,n+1);

    ul=zeros(S*Np,1);

    for s=1:S

        k = 1i*Lambda(s)./c;
        
        
        
        A = 1i/4*besselh(0, k*sqrt((x1-y1).^2 + (x2-y2).^2));
        if min(min(abs(A))) < 1e-20
            A(abs(A) < 1e-20) = 0;
            A = sparse(A);


        end

        ul(idx_Np(s))=A\gl(idx_N(s));    
        
    end
    
    
    phip_hlf(:,n+1) = ul;

hbar.iterate(1);
    
        
end
toc

close(hbar);

phip_hlf(:,M+1-(1:floor((M-1)/2))) = conj(phip_hlf(:,2:floor((M-1)/2)+1));

%% Evaluate the field

Xlim = [-3 3];
Ylim = [-3 3];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*20));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*20));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));



idx=@(s) (s-1)*Nx*Ny+1:(s*Nx*Ny);
idy=@(s) (s-1)*Np+1:(s*Np); 


u_hlf = zeros(S*Nx*Ny,M);


hbar = parfor_progressbar(floor(M/2),['Computing u... M = ',num2str(M/2)]);
tic
for n=0:floor(M/2)

    [P,Lambda]=eig(Am1-lambda*zN^(-(n))*B);
    Lambda=diag(Lambda)/dt;    
    
    gl=phip_hlf(:,n+1);
      
    ul=zeros(S*Nx*Ny,1);
    
    for s=1:S
        
        
    u_s = zeros(Nx*Ny, 1);
    
    
    k1 = 1i*Lambda(s)/c;
    mu = gl(idy(s));
    for m = 1:Np
        u_s = u_s+ mu(m)*1i/4*besselh(0, k1*sqrt((pts(:,1)-xp(m, 1)).^2 + (pts(:,2)-xp(m, 2)).^2));
    end
    
    
    ul(idx(s))=u_s;
    
    end 
    
    
    u_hlf(:,n+1)=kron(P,speye(Nx*Ny))*ul;
    hbar.iterate(1);

    
end
toc
% 
close(hbar);
u_hlf(:,M+1-(1:floor((M-1)/2))) = conj(u_hlf(:,2:floor((M-1)/2)+1));
%%
% Inverting Z-transform
Lam = repmat(lambda.^(0:M-1),S*size(pts,1),1);

u = Lam.^(-1).*ifft(u_hlf,[],2);

u = u((S-1)*Nx*Ny+1:end,:);

u = [zeros(size(u, 1), 1) u];

uinc = zeros(size(u));
for n=1:M

[uinc(:,n),~] = incident_field(pts,c(1)*(tt(n)));

end

u(~in, :) = u(~in, :) - uinc(~in, :);
u(in,:)= 0;


%% Plot the solution
close all  


X = reshape(pts(:,1),Ny,Nx);

Y = reshape(pts(:,2),Ny,Nx);


%%
count = 1;
for n = 1:M
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx), 'EdgeColor', 'black'); hold on;
    plot3([x(:,1); x(1,1)],[x(:,2); x(1,2)],1*ones(N+1,1),'LineWidth',2,'color','k'); 
    plot3([xp(:,1); xp(1,1)],[xp(:,2); xp(1,2)],1*ones(Np+1,1),'LineWidth',2,'color','r'); 
    shading interp;
    hold off
    view(2)
    
    colormap('jet')
    caxis([-1 1])
    axis equal
    axis tight
    axis off
   
    if n < 10
    saveas(gcf, ['fig-beamer/ellipse/ellipse-00',num2str(count),'.png']);
    

    
    elseif n < 100
    saveas(gcf, ['fig-beamer/ellipse/ellipse-0',num2str(count),'.png']);
    
    
    else 
    saveas(gcf, ['fig-beamer/ellipse/ellipse-',num2str(count),'.png']);
        
    end

    drawnow
    
    pause(dt)
    
    count = count +1;

end


