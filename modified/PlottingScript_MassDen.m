<<<<<<< HEAD

% this script opens files containing Mass density data, reads data into 
%   arrays, plots the data in figure, saves the figure, then closes files

Tot_X_Pts = 61;
n = 16;

% open files for reading

fileID_D = fopen( 'MrhoDvals', 'r');

fileID_E = fopen( 'MrhoEvals', 'r');

% size arrays to store data

sizeMrhoD = [Tot_X_Pts,1];

sizeMrhoE = [Tot_X_Pts,1];

% read in data to arrays

MrhoDvals = fscanf(fileID_D, '%f', sizeMrhoD);
MrhoEvals = fscanf(fileID_E, '%f', sizeMrhoE);

% set up x-axis points for plots

x = linspace(0, 3, Tot_X_Pts);

% plot exact results w/o line

plot(x, MrhoEvals, '-o', 'LineStyle', 'none');

hold on;

% plot data

plot(x, MrhoDvals);

% write figure title and label axis

title('Mass density vs. nozzle position');
xlabel('Nozzle position x');
ylabel('Mass density \rho');
legend('Exact', 'Simulation');
%axis([0 3 0.55 0.59])

fclose(fileID_D);
fclose(fileID_E);
=======

% this script opens files containing Mass density data, reads data into 
%   arrays, plots the data in figure, saves the figure, then closes files

Tot_X_Pts = 61;
n = 16;

% open files for reading

fileID_D = fopen('MrhoDvals', 'r');

fileID_E = fopen('MrhoEvals', 'r');

% size arrays to store data

sizeMrhoD = [Tot_X_Pts,1];

sizeMrhoE = [Tot_X_Pts,1];

% read in data to arrays

MrhoDvals = fscanf(fileID_D, '%f', sizeMrhoD);
MrhoEvals = fscanf(fileID_E, '%f', sizeMrhoE);

% set up x-axis points for plots

x = linspace(0, 3, Tot_X_Pts);

% plot exact results w/o line

plot(x, MrhoEvals, '-o', 'LineStyle', 'none');

hold on;

% plot data

plot(x, MrhoDvals);

% write figure title and label axis

title('Mass density vs. nozzle position');
xlabel('Nozzle position x');
ylabel('Mass density \rho');
legend('Exact', 'Simulation');
%axis([0 3 0.55 0.59])

hold off

fclose(fileID_D);
fclose(fileID_E);
>>>>>>> f490d855760e40fead7062a879481d9ddf9d79b7