function [g,grad_g] = incident_field(x,t)
% Planewave incident field
% g(x \cdot d - ct)



a =  -pi/2; % angle of incidence
lag = 4;    % time lag

d = [cos(a) sin(a)];


dor_x_alpha =x(:,1)*d(1)+x(:,2)*d(2); % x \cdot d

beta = 5;
omega = 4.;    % angular frequency
H = @(y) exp(beta*y) ./ (1 + exp(beta*y));  

% Hp = @(y) beta* exp(beta*y) ./ ((1+ exp(beta*y)).^2);

l = 4;

HH = @(y) H(y) .* H(l-y);

% HH = @(y) H(y);
% % HHp = @(y) Hp(y);


g = sin(omega*(t-lag-dor_x_alpha)) .* HH(t-lag-dor_x_alpha);

grad_g = 0;
% grad_g(:,1) = -d(1)*omega*cos(omega*(t-lag-dor_x_alpha)) .* HH(t-lag-dor_x_alpha) ...
%     +sin(omega*(t-lag-dor_x_alpha)) .* d(1).*HHp(t-lag-dor_x_alpha);
% grad_g(:,2) = -d(2)*omega*cos(omega*(t-lag-dor_x_alpha)) .* HH(t-lag-dor_x_alpha) ...
%     +sin(omega*(t-lag-dor_x_alpha)) .* d(2).*HHp(t-lag-dor_x_alpha);
