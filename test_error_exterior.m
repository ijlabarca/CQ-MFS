%% Method of Fundamental Solutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ignacio Labarca 
% Seminar for Applied Mathematics, ETH Zurich
% email:  ignacio.labarca@sam.math.ethz.ch
% date:   July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %     Frequency Domain test

% %     We test the MFS for different complex wavenumbers 'ks'
% %     and number of charge points 'Nps' for a given number of collocation
% %     points 'N', radius R and Rp.
% % 
% %     The test consist in solving the exterior scattering problem at a 
% %     unit disk for an incident field consisting in a point source wave
% %     from the interior of the scatterer.
% % 
% %     The solution is the incident field in the exterior domain 
% %     u(x, t) = uinc(x, t).


clear
clc

R = 1;
Rp = 0.8;
N = 600;
Nps = [10 20 40 80 100 140 180 250 500];
ks = [1i 10 10+10i 10i 100 100+100i 300+300i 500+500i];



theta = -pi;
src = [0.2 0.3]; 


for kss = 1:size(ks, 2)
count = 1;
for Np = Nps
k = ks(kss);


uinc = @(x) 1i/4*besselh(0, k*sqrt((x(:, 1)-src(1)).^2 +(x(:, 2)-src(2)).^2));


t = linspace(0, 2*pi, N).';
tp = linspace(0, 2*pi, Np).';

x = R*[cos(t) sin(t)];
xp = Rp*[cos(tp) sin(tp)];

x1 = repmat(x(:, 1),1, Np);
x2 = repmat(x(:, 2),1, Np);

y1 = repmat(xp(:, 1).',N, 1);
y2 = repmat(xp(:, 2).',N, 1);


%% Least Squares 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 1i/4*besselh(0, k*sqrt((x1-y1).^2 + (x2-y2).^2));
b = -uinc(x);
sol = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Computing total field

% We compute the solution in a grid defined by Xlim and Ylim
Xlim = [-3 3];
Ylim = [-3 3];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*40));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*40));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
in = inpolygon(pts(:, 1), pts(:, 2), x(:, 1), x(:, 2));


u = zeros(Nx*Ny,1);
for n = 1:Np
    u = u+ sol(n)*1i/4*besselh(0, k*sqrt((pts(:,1)-xp(n, 1)).^2 + (pts(:,2)-xp(n, 2)).^2));
end

uinc = 1i/4*besselh(0, k*sqrt((pts(:,1)-src(1)).^2 +(pts(:,2)-src(2)).^2));

u(in) = nan;
uinc(in) = nan;
%% Compute the errors

disp('----------------------------------------')
disp('Relative Error')
disp(count);
disp(norm(u(~in)+uinc(~in))/norm(uinc(~in)))
disp('----------------------------------------')
error(count, kss) = norm(u(~in)+uinc(~in))/norm(uinc(~in));
count = count+1;
end % Nps

end % ks

%% Convergence figure

legend_plot = cell(1, size(ks,2));
for n = 1:size(ks,2)
   legend_plot{n} = "$k = " + string(ks(n)) + "$"; 
end

loglog(Nps, error(:, 1), '--.', 'Linewidth', 2,'MarkerSize', 15);hold on;
for n = 2:size(ks,2)
loglog(Nps, error(:, n), '--.', 'Linewidth', 2, 'MarkerSize', 15);
end
legend(legend_plot, 'Interpreter','latex','FontSize', 14,'Location','SouthWest');
set(gca,'FontSize', 14);
title({'Convergence of the exterior', 'frequency domain problem in a disk'},'Interpreter', 'latex', 'FontSize', 24)
xlabel('Number of charge points','Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative error', 'Interpreter', 'latex', 'FontSize', 18);
yticks([1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1]);
ylim([1e-16 1]);
xlim([10 500]);
xticks(Nps);hold off

    
