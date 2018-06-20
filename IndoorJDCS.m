% Cristina Gómez & James Serna
% Indoor positioning using Joint Distributed Compressed Sensing
% Framework and algorithms based on the proposal by
% "Localization in wireless networks based on jointly compressed sensing"
% by S. Nikitaki, P. Tsakalides. IEEE EUSIPCO, 2011.
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%
% General Framework
%%%%%%%%%%%%%%%%%%%
%  Grid with D cells, K number of mobile users (we start fixing K=1<<<D
%  as in the paper), J Access Points, N time snapshots.

K = 1; % number of targets
J = 10; % number of access points
N = 100; % number of time snapshots
freq = 2.4e3; % frequency in MHz
lambda = 2*pi/freq; % wavelength
M = 10; %number of time samples taken for the CS approach, M/N=compression_rate
snr = -100; %en dB

% Building the grid
D = 100; % number of cells in the grid, must be squared i.e.: 10x10 then D=100.
num_x_y = sqrt(D); % number of columns and rows in the grid, assuming it is squared
area_grid = 100; % area covered of the grid in mts^2
area_cell = round(area_grid/D); %area for each cell
central_cell_0 = [area_cell/2, area_cell/2]; % punto central primer celda, de ahí construimos las demás
cell_center = [];
for xx=0:num_x_y-1
    for yy=0:num_x_y-1
        cell_center = [cell_center; xx+central_cell_0(1) yy+central_cell_0(2)];
    end
end
figure(1)
plot(cell_center(:,1),cell_center(:,2),'kx');
hold on; grid on;

%%%%%%%%%%%%%%%%%%%
% 1. Training Phase
%%%%%%%%%%%%%%%%%%%
% In this phase, each Access Point builds a dictionary of RSSI measurements
% for each position of the grid. We will have then a total of J Psi_j
% matrices, each one of size (NxD)

% Elección aleatoria de posiciones de los AP en el área
% index_position_AP = randi([1 size(cell_center,1)],[J,1]); %falta garantizar que no elija posiciones repetidas,
                                                            %mirar la función ismember para esto

% Por ahora la dejamos fija para pruebas de funcionamiento:
index_position_AP = [1 10 91 100 83 75 32 18 27 63]';

position_AP = cell_center(index_position_AP,:);
figure(1)
plot(position_AP(:,1),position_AP(:,2),'ro');

for jj=1:J
    %OJO! en este cálculo de distancias aún no se tiene contemplado que no
    %puede estar en donde está el AP, por lo cual a veces da una distancia
    %de 0 que luego arroja resultados tipo NaN, pero
    %que necesitaremos revisar si funciona bien!!!
    distances(jj,:) = (sqrt((cell_center(:,1)-position_AP(jj,1)).^2+(cell_center(:,2)-position_AP(jj,2)).^2))/1000; %distance in Kms for RSSI/PL Friis
    Psi_j_one_sample = 30+30-(32.45+20*log10(freq)+20*log10(distances(jj,:)));
    Psi_j_multiple_samples = awgn(repmat(Psi_j_one_sample,N,1),snr);
    Psi_j{jj} = Psi_j_multiple_samples;
    Psi_j{jj}(Psi_j{jj}==Inf) = 0; %esto es para evitar que al comparar con la posición del AP genera luego un NaN
end

%%%%%%%%%%%%%%%%%%
% 2. Runtime Phase
%%%%%%%%%%%%%%%%%%
% In this phase the mobile transmits a beacon that is captured by the APs
% in order to estimate its location

% Elección aleatoria de la posición de los targets en el área
% index_target_position = randi([1 size(cell_center,1)],[K,1]); %falta garantizar que no elija posiciones repetidas,
                                                          %mirar la función ismember para esto
% Por ahora se fija con K=1
index_target_position = 55;

target_position = cell_center(index_target_position,:);

figure(1)
plot(target_position(:,1),target_position(:,2),'b^')

% SIIII FUNCIONA:
% generating the measurements yj=Phi_j*xj=Phi_j*Psi_j*b, partimos de
% construir el b verdadero, este es el que se va a recuperar y lo que nos
% debería dar como resultado el algoritmo de reconstrucción
b = zeros(D,1);
b(index_target_position) = 1; %esta es la posición real del target node
for jj=1:J
    Phi_j{jj} = (1/D)*randn(M,N);   %random Gaussian measurement matrix, from Gaussian Distribution (0,1/D)
    Theta_j{jj} = Phi_j{jj}*Psi_j{jj};
    y_j{jj} = Theta_j{jj}*b;
end

%luego probando cuando yj se genera como una medida RSSI normal, de la cual
%se toman M<<N medidas en el tiempo
% for jj=1:J
%     distance_target(jj,:) = distances(jj,index_target_position);
%     x_j_one_sample = 30+30-(32.45+20*log10(freq)+20*log10(distance_target(jj,:)));
%     x_j{jj} = awgn(repmat(x_j_one_sample,N,1),snr);
%     Phi_j{jj} = (1/D)*randn(M,N);   %random Gaussian measurement matrix, from Gaussian Distribution (0,1/D)
% %     Phi_j{jj} = randi([0 1],M,N);  %trying Bernoulli
%     Theta_j{jj} = Phi_j{jj}*Psi_j{jj};
%     y_j{jj} = Phi_j{jj}*x_j{jj};
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALGORITMO DE RECONSTRUCCIÓN: DCS-SOMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_prima = zeros(D,1);
b_prima = zeros(D,1); %inicializamos en ceros

for jj=1:J
    b_j{jj,:} = Theta_j{jj}'*y_j{jj};
    b_prima = b_prima + abs(b_j{jj});    
end
[max_b,p_l] = max(b_prima);
plot(cell_center(p_l,1),cell_center(p_l,2),'g*')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGUNOS ANÁLISIS DE SENSIBILIDAD HASTA AHORA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Con SNR muy alta se queda estático pero en valores errados, poner valores
% similares a los de referencia que ellos dan en el artículo, -70, -90... y
% empieza a funcionar
% La topología de la red, el número de AP y su ubicación afectan el
% comportamiento
% La tasa de compresión también afecta que funcione o no el algoritmo
