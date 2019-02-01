% ==============================================================================
% Function to calculate velocity fields for a 2 dimensional field
% (it is external so different velocity fields can be tested
% easily) Antonio Preziosi-Ribero, November 2017
% Universidad Nacional de Colombia - Northwestern University
% ==============================================================================

function [u v] = velocity(x, y, db, uu, vs, qin)

        % Functions described in paper by Packman et al. 2000
        % adding inflow or outflow conditions from experiments

	% With underflow
	% u  = -cos(x) .* (tanh(db) * sinh(y) + cosh(y)) + uu;

	% No underflow
	u  = -cos(x) .* (tanh(db) * sinh(y) + cosh(y));

	% u = -cos(x) .* exp(y);

    u = real(u);

    v = -sin(x) .* (tanh(db) * cosh(y) + sinh(y));

	% v = -sin(x) .* exp(y);

   	v = v + vs + qin;

    v = real(v);

    % Testing functions - linear in x and y direction (1 mm/s
    % in each of the directions)

%        u = zeros(size(x)) * 1;

%        v = ones(size(x)) * -1;
