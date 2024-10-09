% main.m

function main()
	addpath("..");

	% Initial guess for (xB, yB, zB)
	xB_initial = 1; % meters
	yB_initial = 1; % meters
	zB_initial = 1; % meters
	tol = 1e-8; % tolerance


	initial_guess = [xB_initial; yB_initial; zB_initial]; % vector of xB, yB and zB

	% Equation 1, Equation 2 and Equation 3 | v(1) = xB, v(2) = yB v(3) = zB
	Eq{1} = @(v) v(1)^2 - 2*v(2)^2 - v(2) - 2*v(3);
	Eq{2} = @(v) v(1)^2 - 8*v(2)^2 + 10*v(3);
	Eq{3} = @(v) (v(1)^2)/(7*v(2)*v(3))-1;

	% Call Newton's method (fn_newthon_method_n3.m)
	[xB, yB, zB] = fn_newton_method_n3(Eq, initial_guess, tol);

	printf("\n==============\nxB: %f\nyB: %f\nzB: %f\n", xB, yB, zB);
end

main();
