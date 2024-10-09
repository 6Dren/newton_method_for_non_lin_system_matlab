% main.m

function main()
	addpath("..");

	% Initial guess for (xB, yB, zB)
	xB_initial = 0.1; % meters
	yB_initial = 0.1; % meters
	zB_initial = -0.1; % meters
	tol = 1e-8; % tolerance


	initial_guess = [xB_initial; yB_initial; zB_initial]; % vector of xB, yB and zB

	% Equation 1, Equation 2 and Equation 3 | v(1) = xB, v(2) = yB v(3) = zB
	Eq{1} = @(v) 3*v(1) - cos(v(2)*v(3)) - 1/2 ;
	Eq{2} = @(v) v(1)^2 - 81*(v(2)^2 + 0.1)^2 + sin(v(3)) + 1.06;
	Eq{3} = @(v) -v(1) + e^(-v(1)*v(2)) + 20*v(3) + (10*pi - 3)/3;

	% Call Newton's method (fn_newthon_method_n3.m)
	[xB, yB, zB] = fn_newton_method_n3(Eq, initial_guess, tol);

	printf("\n==============\nxB: %f\nyB: %f\nzB: %f\n", xB, yB, zB);
end

main();
