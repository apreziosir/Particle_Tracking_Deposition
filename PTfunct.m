% ==============================================================================
% PARTICLE TRACKING FUNCTION TO ESTIMATE REMOBILIZATION OF PARTICLES
% Antonio Preziosi Ribero (July 2018)
% Fine particles' deposition in a single bedform modeled with losing and 
% gaining conditions to be compared with experimental data from Shai Arnon and 
% Aryeh Fox (Fox et al. 2014 & 2018 WRR 25014)
% Non-dimensional form of the equations implemented
% Code running in octave with Matlab exceptions in comments
% Universidad Nacional de Colombia - Northwestern University
% ==============================================================================

% close all
% clear all

function GWFL = PTfunct(qin, k0)

	% ==========================================================================
	% Definition physical constants (cgs)
	% ==========================================================================

	rho = 1.0;				% Water density [g/cm^3]
	g = 981.0;				% Gravity [cm/s^2]
	mu = 0.1;				% Dyn viscosity [cm^2/s]

	% ==========================================================================
	% Data from lab - gathered from Arnon & Fox experiments 2017
	% Data not reported found in Boano & Fox paper
	% ==========================================================================

	H = 1.5;				% Dune height [cm]
	d = 8.75;				% Mean stream depth [cm]
	db = 20.0;	  			% Depth of impervious layer [cm]
	U = 17.143;	  			% Mean stream velocity [cm/s]
	Wl = 15.0;	  			% Dune wavelength [cm]
	% qin = 0.00;			% w/l discharge [cm/s]
	phi = 0.33;	  			% Bed porosity [nondim]
	K = 0.1195;	  			% Hyd conductivity [cm/s]
	D50 = 0.0384;			% Sand grain size [cm]
	vs = 0.00;	  			% Filt. velocity [cm/min]
	Area = 1.566;			% Bed area [m^2]
	S = 0.000147;			% Slope [m/m]

	% ==========================================================================
	% Numerical model parameters
	% ==========================================================================

	n = 2e5;				% Number of particles
	T = 1;          		% Total time simulated [-]
	nT = 4e5;   			% Number of timesteps [-]
	sts = 100;				% Save every sts timesteps [-]
	offset = 1e-4;	      	% Part. offset boundary [-]
	poff = 0.55;			% Max portion of the domain filled
	cellsz = 1e-3; 			% Cell size for interpolating [-]

	% Maximum specific uptake rate for filtering  (cm-1)
	% Typical value taken from Packman et. al 2001

	% k0 = 0.6;        			% Filtering coeff. [1/cm]

	% Tuning parameters for the num model

	% interp methods: 'linear', 'nearest', 'spline' (for interp2)
	% No need to use itp_m with qinterp2

	itp_m = 2;				% Interp. method (deprecated)
	tuneKbed = 0.0;				% Tuning parameter in bed

	% ==========================================================================
	% Calculation of parameters for the model based on lab data -  This part is
	% moved up in order to estimate the non dimensional quantities of the model
	% ==========================================================================

	% Calculation of mean head over bedform (outside function)
	hm = h_m(U, H, d, g);

	% Calculation of k (term inside trigonometric functions)
	k = (2 * pi) / Wl;

	% Maximum pumping velocity
	um = k * K * hm;

	% Converting settling  vel and affecting by porosity
	% Modified january 2018 (affecting by porosity) as in Packman's paper
	vs = (vs / 60) * phi;
	disp('Settling velocity (cm/s)[Darcy velocity]:');
	disp(vs);

	% Converting and estimating w/l flow
	% qin = Qin(qin, phi);
	% disp('Darcy inflow/outflow vel (cm/s);')
	% disp(qin)

	% ==========================================================================
	% Transforming to non dimensional quantities - Added June 2018
	% ==========================================================================

	% Non dimensional velocity for estimated lab quantities

	% qin = qin / um; 				% In/out nondimensional velocity
	uu = K * S / um;				% Underflow nondimensional velocity
	vs = vs / um;					% Settling non dimensional velocity
	% k0 = k0 * Wl / (2 * pi) 		% Non dimensional filt. coefficient [-]

	% Non dimensional total time to simulate
	T = k ^ 2 * K * hm * (T * 3600) / phi;

	% Calculating the number of timesteps needed
	dT = T / nT;				% nondimensional timestep size [-]

	% Maximum y (in non dimensional form). It is also the db* for the non
	% dimensional form from Packman 2000.
	ymax = 2 * pi * db / Wl;

	% Refinement for interp2
	refx = floor(2 * pi / cellsz); 		% x refinement (portion of x length)
	refy = floor(ymax / cellsz); 		% y refinement (portion of x length)

	% Generating mapping vectors and meshgrid
	x = linspace(0, 2 * pi, refx);		% x vector of coords [-]
	y = linspace(-ymax, 0, refy);		% y vector of coords [-]
	[X Y] = meshgrid(x, y);				% Meshgrid generation for calculations

	% Initial position of particles - The Y position is affected by machine 
	% precision for interpolating purposes (if particles are seeded in the 
	% boundary the interpolation algorithm doesn't work at all)
	% offst = 2 * pi * offset;
	posx = linspace(0 + offset , pi - offset, n); 	% Initial x pos. [-]
	posy = ones(1,n) * (0 - 100000 * eps);    		% Initial y pos. [-]

	% ==========================================================================
	% Setting up the domain and seeding the particles equally spaced. These are 
	% the initial conditions of the model
	% ==========================================================================

	% Vectors that allocate filtered particles - start with NaN since 0
	% values are plotted.
	posxfil = nan(1, n);				% Filtered x pos.[-]
	posyfil = nan(1, n);				% Filtered y pos.[-]
	nfil = 0;			        		% Number of filt particles

	% Vector that counts amount of remobilized particles in one
	% timestep
	remob = zeros(1, nT);

	% ==========================================================================
	% Parameters for saving information, this setup includes data and
	% video parameters.
	% ==========================================================================

	% Save information every sts timesteps
	disp('Saving information every [timesteps]:');
	disp(sts);

	% Setting up arrays that are going to save information and saving
	% initial information
	keepind = 1;
	keepframes = 0:sts:nT;
	keepx = nan(n, length(keepframes));
	keepy = nan(n, length(keepframes));
	keepxfil = nan(n, length(keepframes));
	keepyfil = nan(n, length(keepframes));
	keepx(:,1) = posx;
	keepy(:,1) = posy;

	% ==========================================================================
	% Calculation of velocity fields u and v with x and y vectors (the
	% refinement of the calculation will affect the interpolation in the
	% loop). Also, calculation of the Hydraulic conductivity fields
	% Klong and Kvert.
	% ==========================================================================

	% Estimating the value of the velocity fields (outside functions)
	[u v] = velocity(X, Y, ymax, uu, vs, qin);

	% Estimating the conductivity fields value
	[Klong Kvert] = hyd_cond(X, Y, phi, D50, K);

	% ==========================================================================
	% Temporal loop of the program
	% ==========================================================================

	tic
	for t = 1:nT

		% Obtaining interpolated velocities for each position qinterp2
		uinterp = qinterp2(X, Y, u, posx, posy, itp_m);
		vinterp = qinterp2(X, Y, v, posx, posy, itp_m);

		% Obtaining interpolated hydraulic conductivity for each pos qinterp2
		if tuneKbed ~= 0
		    kvinterp = qinterp2(X, Y, Kvert, posx, posy, itp_m);
		    klinterp = qinterp2(X, Y, Klong, posx, posy, itp_m);
		end

		% New x coordinate for timestep t
		if tuneKbed ~= 0

		    tempposx = posx + (uinterp * dT) + (((2 * klinterp * dT) .^ 0.5) ...
		    .* randn(1, n)) * tuneKbed;

		else

		    tempposx = posx + (uinterp * dT);

		end

		% Applying periodic boundary condition as the modulo of the x
		% length of the domain
		tempposx = mod(tempposx, 2 * pi);

		% New y coordinate for timestep t
		if tuneKbed ~= 0

		    tempposy = posy + (vinterp * dT) + (((2 * kvinterp * dT) .^ 0.5) ...
		    .* randn(1, n)) * tuneKbed;

		else

		    tempposy = posy + (vinterp * dT);

		end

		% calculating ds
		ds = dT .* sqrt(uinterp .^ 2 + vinterp .^ 2);

		% Updating posx and posy (pos0 = pos1)
		posx = tempposx;
		posy = tempposy;

		% Applying vertical boundary conditions (bounce inside the
		% domain) Same as Angang code
		low = find(posy < -ymax);
		posy(low) = -ymax - posy(low);

		% Applying top BC (remobilized particles)
		high = find(posy > 0);
		remob(t) = length(find(posy > 0));

		% Taking particles out of the domain
		posy(high) = NaN;
		posx(high) = NaN;

		% Specify filtering part of the model - Selecting particles
		% that are going to be filtered
		k0interp = ones(1,n) * k0;
		k0interp(posy < 0) = k0;
		filtered  = rand(1, n) - (k0interp .* ds);

		% Saving the filtered positions in separate arrays
		posxfil(filtered <= 0) = posx(filtered <= 0);
		posyfil(filtered <= 0) = posy(filtered <= 0);

		% Taking out particles from calculations
		posx(filtered <= 0) = NaN;
		posy(filtered <= 0) = NaN;

		% Saving positions in selected timesteps
		if ismember(t, keepframes)

			filtradas = length(posx(posx == NaN));
			% Displaying relevant information for every moment that the program
			% saves information
	%		disp('==========================================')
	%	    disp(['Saving data corresponding to time (-):']);
	%	    disp(t * dT);
	%		disp('Remobilized particles:');
	%		disp(remob(t));
	%		disp('Filtered particles at this moment:')
	%		disp(filtradas)
	%		disp('==========================================')

			% Saving information to vectors for every snapshot
		    keepind = keepind + 1;
		    keepx(:,keepind) = posx;
		    keepy(:,keepind) = posy;
		    keepxfil(:,keepind) = posxfil;
		    keepyfil(:,keepind) = posyfil;

		end

		% Total time simulated until now
		tot_t = t * dT;

		% Checking for particle filtration and remobilization (if all particles
		% are NaN, then the loop can stop here)
		exit_loop = length(find(posx == NaN));

		if exit_loop == n

			disp('All particles have been filtered or remobilized');
			break
		end

	% End of the time loop (particle tracking finished)
	end

	% Displaying time consumed for particle tracking (checking in screen)
	disp('Time consumed for estimating particle tracking');
	toc

	GWFL.Lx = 2 * pi;
	GWFL.Ly = ymax;
	GWFL.n = n;
	GWFL.sts = sts;
	GWFL.dT = dT;
	GWFL.X = X;
	GWFL.Y = Y;
	GWFL.u = u;
	GWFL.v = v;
	GWFL.keepframes = keepframes;
	GWFL.keepx = keepx;
	GWFL.keepy = keepy;
	GWFL.keepxfil = keepxfil;
	GWFL.keepyfil = keepyfil;
	GWFL.remob = remob;
