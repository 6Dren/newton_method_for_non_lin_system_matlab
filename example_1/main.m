% main.m

function main()
	% Initial guess for (xB, yB)
	xB_initial = 1500; % meters
	yB_initial = 2000; % meters
	tol = 1e-8; % tolerance

	% m (Note: dnt confuse meters(m) with kilometers(1 km = 1e3 = 10³ m))
	d1 = 3e3; % meters
	d2 = 4e3; % meters

	% seconds
	t1 = 0; 
	t2 = 5.72e-6;
	t3 = 8.58e-6;

	initial_guess = [xB_initial; yB_initial]; % vector of xB and yB

	Eq = equations(d1, d2, t1, t2, t3);

	% Call Newton's method (fn_newthon_method_non-lin-system.m)
	[xB, yB] = fn_newton_method_non_lin_system(Eq, initial_guess, tol);

	printf("xB = %f, yB = %f\n", xB, yB);
end


function Eq = equations(d1, d2, t1, t2, t3)
	c = 3e8; % speed of light (m/s)

	% Equation 1 and Equation 2 | v(1) = xB and v(2) = yB
	Eq{1} = @(v) (4*(v(1) - d1/2)^2)/(c^2*(t2-t1)^2) - (4*(v(2)^2))/(d1^2 - c^2*(t2-t1)^2) - 1;
	Eq{2} = @(v) (4*(v(2) - d2/2)^2)/(c^2*(t3-t2)^2) - (4*(v(1) - d1)^2)/(d2^2 - c^2*(t3-t2)^2) - 1;
end

main();
