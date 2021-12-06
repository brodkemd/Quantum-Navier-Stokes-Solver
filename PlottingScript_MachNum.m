
% this script opens files containing Mach number data, reads data into 
%   arrays, plots the data in figure, saves the figure, then closes files

Tot_X_Pts = 61;
n = 16;

% open files for reading

fileID_D = fopen('MachDvals', 'r');

fileID_E = fopen('MachEvals', 'r');

% size arrays to store data

sizeMachD = [Tot_X_Pts,1];

sizeMachE = [Tot_X_Pts,1];

% read in data to arrays

MachDvals = fscanf(fileID_D, '%f', sizeMachD);
MachEvals = fscanf(fileID_E, '%f', sizeMachE);

% set up x-axis points for plots

x = linspace(0, 3, Tot_X_Pts);

% plot exact results w/o line

plot(x, MachEvals, '-o', 'LineStyle', 'none');

hold on;

% plot data

plot(x, MachDvals);

% write figure title and label axis

title('Mach number vs. nozzle position');
xlabel('Nozzle position x');
ylabel('Mach number M');
legend('Exact', 'Simulation');
%axis([0 3 0.55 0.59])

fclose(fileID_D);
fclose(fileID_E);
