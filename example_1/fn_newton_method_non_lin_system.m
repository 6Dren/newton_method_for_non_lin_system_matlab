% fn_newton_method_non-lin-system.m
% New's Method for system of Non-linear Equations (fn_newton_method_non-lin-system.m)

% In this example can be use in the LORAN's Formula
% Eq 1 | (4(xB - d1/2)²)/(c²(t2 - t1)²) - (4yB²)/(d1² - c²(t2 - t1)²) = 1
% Eq 2 | (4(yB - d2/2)²)/(c²(t3 - t2)²) - (4(xB - d1)²)/(d2² - c²(t3 - t2)²) = 1

% where:
% eq 1, eq 2: Equation 1 and Equation 2 (bcs it's a system)
%  xB, yB are the position that we looking for
%  d1, d2 are the distances between three points
%  t1, t2, t3 are the time's signal
%  c is the speed of light: 3 * 10⁸ (approximately)

function [x, y] = fn_newton_method_non_lin_system(Eq, initial_guess, tol)
	c = 0; % counter
	pos = initial_guess; % position

	err = tol + 1;

	while err > tol
		J = numerical_jacobian(Eq, pos); % J: jabobian of the equation (Eq(1) and Eq(2)) with x starting with the initial guess
		
		Eq_val = [Eq{1}(pos); Eq{2}(pos)];

		if abs(det(J)) < 1e-10 % the determinant can't be zero!
			error("Jacobian is singular or nearly singular!\n");
		end

		% update position
		delta = -J \ Eq_val; % J⋅δ=−F => δ = -J⋅\ F     ' \ ' solves system of equations 
		pos = pos + delta;
		err = norm(delta);

		% Show to the user a counter of the loop with all variables used
		printf("\n\n---------------------\nCounter: %d\n", ++c);
		printf("---------------------\n");
		disp("J: ");
		disp(J);
		printf("---------------------\n");
		disp("Eq_val: ");
		disp(Eq_val);
		printf("---------------------\n");
		disp("| det(J) |: ");
        disp(abs(det(J)));
		printf("---------------------\n");
		disp("delta: ");
        disp(delta);
		printf("---------------------\n");
		disp("position (x, y): ");
        disp(pos);
		printf("---------------------\n");
		disp("error: ");
        disp(err);
	end

	printf("\n");
	x = pos(1);
	y = pos(2);
end

function J = numerical_jacobian(Eq, pos)
	% Numerical differentiation
	h = 1e-6; % since we can't divide by zero, we choice a value pretty close to 0
	J = zeros(2,2);

	% Calculate the Equations at the current position
	Eq1 = Eq{1}(pos); % equation 1
	Eq2 = Eq{2}(pos); % equation 2

	%% Calculate the partial derivatives numerically
	% 1st Equation
	J(1, 1) = (Eq{1}([pos(1) + h, pos(2)]) - Eq1) / h;
	J(1, 2) = (Eq{1}([pos(1), pos(2) + h]) - Eq1) / h;

	% 2st Equation
	J(2, 1) = (Eq{2}([pos(1) + h, pos(2)]) - Eq2) / h;
	J(2, 2) = (Eq{2}([pos(1), pos(2) + h]) - Eq2) / h;
end
