% ==============================================================================
% Function to estimate the mean value of head based on lab
% data. Taken from Packman et al. 2000 WRR (which is taken from
% Elliott and Brooks 1997). 
% Antonio Preziosi-Ribero, November 2017
% Universidad Nacional de Colombia - Northwestern University
% ==============================================================================

function h = h_m(U, H, d, g)

	if H / d <= 0.34
            
			h = 0.28 * (U^2 / (2 * g)) * (H / (d * 0.34)) ^ (3/8);
            
        else
            
			h = 0.28 * (U^2 / (2 * g)) * (H / (d * 0.34)) ^ (3/2);
            
end
