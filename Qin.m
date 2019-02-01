% ==============================================================================
% Function to estimate in/out velocity from data gathered in the
% lab. This is an enhancing made to simplify the code that is
% already running.
% This function is not implemented in the code but it is provided if needed
% Antonio Preziosi-Ribero, January 2018
% Universidad Nacional de Colombia - Northwestern University
% ==============================================================================

function qin = Qin(qin, phi)

	% qin in cm/s - from cm/d to cm/s
	qin = qin * (1 / 86400) * (1 / phi);

    
