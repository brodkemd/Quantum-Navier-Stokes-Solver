
% this script opens files containing temperature data, reads data into 
%   arrays, plots the data in figure, saves the figure, then closes files

Tot_X_Pts = 61;
n = 16;

% open files for reading

fileID_D = fopen( ...
'D:\OneDrive - University of Cincinnati 1\research\software\quantum_NS_solver\Gaitan_quantum_Navier_Stokes_simulation_software\TempDvals', 'r');

fileID_E = fopen( ...
'D:\OneDrive - University of Cincinnati 1\research\software\quantum_NS_solver\Gaitan_quantum_Navier_Stokes_simulation_software\TempEvals', 'r');

% size arrays to store data

sizeTempD = [Tot_X_Pts,1];

sizeTempE = [Tot_X_Pts,1];

% read in data to arrays

TempDvals = fscanf(fileID_D, '%f', sizeTempD);
TempEvals = fscanf(fileID_E, '%f', sizeTempE);

% set up x-axis points for plots

x = linspace(0, 3, Tot_X_Pts);

% plot exact results w/o line

plot(x, TempEvals, '-o', 'LineStyle', 'none');

hold on;

% plot data

plot(x, TempDvals);

% write figure title and label axis

title('Temperature vs. nozzle position');
xlabel('Nozzle position x');
ylabel('Temperature T');
legend('Exact', 'Simulation');
%axis([0 3 0.55 0.59])

fclose(fileID_D);
fclose(fileID_E);
