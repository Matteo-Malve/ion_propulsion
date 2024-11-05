%%
clear; close all; clc;

r2 = @(x,y) x.^2+y.^2;
uex = @(r2) exp(-r2./0.0004);

x=linspace(-0.004,0.004,201);
y=linspace(0.,0.004,101);
rr = linspace(0,sqrt(2.)*0.004,101);
rrr = linspace(-5,10,1001);

% Create a meshgrid
[XX, YY] = meshgrid(x, y);

% Compute the function values
R2 = r2(XX, YY);
U = uex(R2);

% Plot the result
subplot(1,3,1)
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
title('Color Plot of uex = exp(-r^2)');
xlabel('x');
ylabel('y');

subplot(1,3,2)

plot(rrr,uex(rrr));
hold on
plot(rr,uex(rr),"Color","red","LineWidth",3);
xline(0)
yline(0)

subplot(1,3,3)

plot(rr,uex(rr));

%%
clear; close all; clc;

uex = @(x,y) 1 - ((x./0.004).^2) .* ((y./0.004));

x=linspace(-0.004,0.004,201);
y=linspace(0.,0.004,101);

% Create a meshgrid
[XX, YY] = meshgrid(x, y);

% Compute the function values
U = uex(XX,YY);

% Plot the result
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
title('Color Plot of uex = exp(-r^2)');
xlabel('x');
ylabel('y');

%%
clear; close all; clc;

Ve = 20000.;
RE = 0.0004;
R0 = 3 * RE;

uex = @(x,y) Ve .* (1 - (x.^2)./(R0.^2) - (y.^2)./(R0.^2)) .* ((x.^2)./(RE.^2) + (y.^2)./(RE.^2) - 1);

% Create a meshgrid
x=linspace(-0.004,0.004,201);
y=linspace(0.,0.004,101);
[XX, YY] = meshgrid(x, y);


% Compute the function values
U = uex(XX,YY);

% Plot the result
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
xlabel('x');
ylabel('y');

%%
clear; close all; clc;

Ve = 20000.;
l = 0.0004;
R = sqrt(2)*l;
sigma2 = 0.01;

uex = @(x,y) exp(-(sqrt(x.^2+y.^2)-R).^2./(2*sigma2));

% Create a meshgrid
x=linspace(-0.004,0.004,201);
y=linspace(-0.004,0.004,201);
[XX, YY] = meshgrid(x, y);


% Compute the function values
U = uex(XX,YY);
mask = sqrt(XX.^2 + YY.^2) < R;
U(mask) = NaN;

% Plot the result
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
xlabel('x');
ylabel('y');



%%
clear; close all; clc;

Ve = 20000.;
l = 0.0004;
R = sqrt(2)*l;
sigma2 = 0.0000005;

uex = @(x,y) exp(-(sqrt(x.^2+y.^2)-R).^2./(2*sigma2)) ...
     .*  Ve .* 0.5 ...
     .* (1 - sin(x./l.*pi) .* sin(y./l.*pi));

% Create a meshgrid
x=linspace(-0.004,0.004,201);
y=linspace(-0.004,0.004,201);
[XX, YY] = meshgrid(x, y);

U = uex(XX,YY);

% Plot the result
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
xlabel('x');
ylabel('y');

hold on;
%plot([-l, l, l, -l, -l], [l, l, -l, -l, l], 'k-', 'LineWidth', 2); % Square outline in black
hold off;
%%
close all;
ray = linspace(l,0.004,401);
same_y = l./2.*ones(size(ray));
plot(ray,uex(ray,same_y))

%%
clear; close all; clc;

Ve = 20000.;
l = 0.0004;
R = sqrt(2)*l;
sigma2 = 0.0000005;

uex = @(r) - Ve .* (r-l) .* (r-0.004);

rr = linspace(l,0.004,401);
plot(rr,uex(rr))

%%
clear; close all; clc;

Ve = 20000.;
l = 0.0004;
L = 0.0004;
R = sqrt(2)*l;
sigma2 = 0.0000005;
eps_0 = 8.854;
eps_r = 1.0006;

uex = @(x,y) - eps_0 * eps_r * Ve .* ( l + L - 4.*sqrt(x.^2+y.^2) ) ./ sqrt(x.^2+y.^2);

% Create a meshgrid
x=linspace(-0.004,0.004,201);
y=linspace(-0.004,0.004,201);
[XX, YY] = meshgrid(x, y);

U = uex(XX,YY);

% Plot the result
imagesc(x, y, U);
colorbar;
axis xy; % Corrects axis orientation
xlabel('x');
ylabel('y');

hold on;
%plot([-l, l, l, -l, -l], [l, l, -l, -l, l], 'k-', 'LineWidth', 2); % Square outline in black
hold off;
