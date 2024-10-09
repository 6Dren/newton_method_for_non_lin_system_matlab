% main.m

function main()
	addpath("..");

	% Initial guess for (xB, yB)
	xB_initial = 1500; % meters
	yB_initial = 2000; % meters
	tol = 1e-8; % tolerance


	initial_guess = [xB_initial; yB_initial]; % vector of xB, yB

	% Equation 1 and Equation 2 | v(1) = xB, v(2) = yB
	Eq{1} = @(v) (v(1)^2)/(186^2) - (v(2)^2)/(300^2 - 186^2) - 1;
	Eq{2} = @(v) ((v(2) - 500)^2)/(279^2) - ((v(1) - 300)^2)/(500^2 - 279^2) - 1;

	% Call Newton's method (fn_newthon_method_n2.m)
	[xB, yB] = fn_newton_method_n2(Eq, initial_guess, tol);

	printf("\n==============\nxB: %f\nyB: %f\n", xB, yB);
end

main();
