
% ==============================================================================
% Code that calls the Particle Tracking function with a given vector of inflow
% outflow conditions and a vector of filtration coefficients.

% Antonio Preziosi-Ribero
% Universidad Nacional de Colombia - Sede Bogotá
% Facultad de Ingeniería
% ==============================================================================

% ==============================================================================
% Edits: 
% 07/2018
% 02/01/2019 
% ==============================================================================

% close all 
% clear all 

% ==============================================================================

% Vector with dimensionless in/out flows (fraction of maximum pumping velocity)
q_io = [-0.08, 0.0, 0.08];

% Vector with dimensionless filtering coefficient 
k0 = [1.45];

% Nested loops to run all the possible conditions imposed
for i = 1:length(k0)

	for j = 1:length(q_io)

		%  Getting info from the particle tracking model
		DATA_1 = PTfunct(q_io(j), k0(i));

		% Changing the name of the mat file to be saved
		tit = strcat('GWFL_', num2str(q_io(j)), '_', num2str(k0(i)), '.mat');

		% Writing mat file to local folder (mat files are readable as hdf5 files
		% in any platform like hdfcompass)
		save(tit, 'DATA_1','-v7.3');

		% Displaying on-screen information
		disp('Termine estos datos:');
		disp(q_io(j));
		disp(k0(i));

	end

end

exit
