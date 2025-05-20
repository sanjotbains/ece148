% Updated MATLAB script with custom axis ticks/labels

% Parameters
omega0 = 2 * pi;               % Sampling frequency
Delta_t = 2 * pi / omega0;     % Sample spacing
omega_x = linspace(-4*omega0, 4*omega0, 1000);  % CTFT domain

% --- Figure 1: CTFT mapping with ω₀ scaling ---
omega_y = mod(omega_x + omega0/2, omega0) - omega0/2;

figure (1);
plot(omega_x/omega0, omega_y/omega0, 'b', 'LineWidth', 1.5);
xlabel('\omega_x (units of \omega_0)');
ylabel('\omega_y (units of \omega_0)');
title('CTFT Mapping: \omega_x to \omega_y');
grid on;
xticks(-4:1:4);
yticks(-0.5:0.25:0.5);
xlim([-4, 4]);
ylim([-0.5, 0.5]);

% --- Figure 2: DTFT mapping with π scaling ---
theta_y = mod(omega_x * Delta_t + pi, 2*pi) - pi;

figure (2);
plot(omega_x/pi, theta_y/pi, 'r', 'LineWidth', 1.5);
xlabel('\omega_x (units of \pi)');
ylabel('\theta_y (\units of\pi)');
title('DTFT Mapping: \omega_x to \theta_y');
grid on;
xticks(-4:1:4);
yticks(-1:0.5:1);
xlim([-4, 4]);
ylim([-1, 1]);
