% ==============================================================================
% Function to estimate the hydrualic conductivity field. In the
% first example this value will be constant. Later the value can be
% modified accordng to filtering results.
% Antonio Preziosi-Ribero, November 2017.
% Universidad Nacional de Colombia - Northwestern University
% ==============================================================================

function [Klong Kvert] = hyd_cond(X, Y, phi, D50, K);

	% Constant conductivity matrix
	Klong = ones(size(X)) * K;
	Kvert = Klong;
    
end