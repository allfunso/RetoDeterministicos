tic
clear
clc
close all

N = 2^9;
L = 30;
dx = L/N;     dk = 2*pi/L;

x = [-N/2:1:N/2-1]*dx;    k = [-N/2:1:N/2-1]*dk;

kshift = fftshift(k);     kshift2 = kshift.^2;

dz = dx.^2/4;
zfinal = 27;
steps = zfinal/dz;
h = dz*((1:1:steps)-1);

steps = ceil(zfinal/dz);

% Límite de A es 3
A = 2.9;
A1 = A;   W1 = 1/A1;
A2 = A;   W2 = 1/A2;
A3 = A;   W3 = 1/A3;
A4 = A;   W4 = 1/A4;
A5 = A;   W5 = 1/A5;

% Límite de q es 5
q = 3.81;
uo = A1*sech((x-2*q)/W1) + A2*sech((x-1*q)/W2) + A3*sech((x)/W3) + A4*sech((x+1*q)/W4) + A5*sech((x+2*q)/W5);
un = uo;
InitialMax = max(abs(uo));
interval = 20; % Cada cuantas iteraciones se guarda un valor y se calcula potencia

% Generar vectores de u y de potencias
utotal = zeros(N, floor(steps/interval));
potencias = zeros(floor(steps/interval), 1);
potencias(1) = sum(uo.*conj(uo) * dx);

for count = 1:1:steps
    F_NL = fft(exp(1i*dz*abs(un).^2).*un);
    F_D = exp(-1i*kshift2*dz/2).*F_NL;
    un = ifft(F_D);
    if mod(count, interval) == 0
        idx = count/interval;
        utotal(:, idx) = un;
        potencias(idx+1) = sum(un.*conj(un) * dx);
    end
end

B = abs(utotal);
B_size = size(B);
v = B_size(2);
idx1 = round(v/3); % Indice xi 1/3
idx2 = round(2*v/3); % Indice x 2/3

figure(1)
imagesc(h, x, B);
title('Propagaci\''on de Pulso Solitario', 'Interpreter','latex', 'FontSize',18);
colormap(hot(256))
c = colorbar;
set(c,'TickLabelInterpreter','latex', 'FontSize',13);

xline(zfinal/3,'w');
xline(2*zfinal/3,'w');

xlabel('Distancia','Interpreter','latex', 'FontSize',15);
ylabel('Tiempo', 'Interpreter','latex', 'FontSize',15)
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlim([0 27]);

figure(2);
plot(x,abs(uo));
hold on
plot(x,abs(un),'r');
%plot(x, B(:, idx1), "r");
%plot(x, B(:, idx2), "r");
str = 'Pulso final';
plot(x, A*ones(N) + 0.05*A, 'k--')
plot(x, A*ones(N) - 0.05*A, 'k--')
legend('Pulso inicial', str, 'Interpreter','latex', 'FontSize',10);
ylabel('Amplitud','Interpreter','latex', 'FontSize',10); 
xlabel('tiempo','Interpreter','latex', 'FontSize',10);
hold off

maximos1 = B(islocalmax(B(:, idx1)), idx1);
maximos1 = maximos1(maximos1 > 0.9*A);
maximos2 = B(islocalmax(B(:, idx2)), idx2);
maximos2 = maximos2(maximos2 > 0.9*A);
maximos3 = B(islocalmax(B(:, end)), end);
maximos3 = maximos3(maximos3 > 0.9*A);

% Errores porcentuales
errores1 = (1 - maximos1 / A) * 100;
errores2 = (1 - maximos2 / A) * 100;
errores3 = (1 - maximos3 / A) * 100;

fprintf("Errores 1: %.4f %% \n", errores1)
fprintf("Errores 2: %.4f %% \n", errores2)
fprintf("Errores 3: %.4f %% \n", errores3)

%% Conservación de la Potencia

% Graficar potencias y líneas de error +- 0.1%
figure(3)
plot(h(1:interval:end), potencias)
hold on
p_inicial = potencias(1);
plot(h(1:interval:end), p_inicial*ones(length(potencias)) + 0.0001*p_inicial, 'k--')
plot(h(1:interval:end), p_inicial*ones(length(potencias)) - 0.0001*p_inicial, 'k--')
title("Conservación de la potencia")
xlabel("xi")
ylabel("Potencia")
legend("Potencia", "error (0.01%)")
xlim([0, 27])

toc