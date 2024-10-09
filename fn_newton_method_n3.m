% fn_newton_method_non-lin-system.m
% New's Method for system of Non-linear Equations (fn_newton_method_non-lin-system.m)

% where:
% eq 1, eq 2 and eq 3: Equation 1, Equation 2 and Equation 3
%  xB, yB, ZB are the positions that we looking for


function [x, y, z] = fn_newton_method_n3(Eq, initial_guess, tol)
	c = 0; % counter
	pos = initial_guess; % position

	err = tol + 1;

	while err > tol
		J = numerical_jacobian(Eq, pos); % J: jabobian of the equation (Eq(1) and Eq(2)) with x starting with the initial guess
		
		Eq_val = [Eq{1}(pos); Eq{2}(pos); Eq{3}(pos)];

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
		disp("position (x, y, z): ");
        disp(pos);
		printf("---------------------\n");
		disp("error: ");
        disp(err);
	end

	printf("\n");
	x = pos(1);
	y = pos(2);
	z = pos(3);
end

function J = numerical_jacobian(Eq, pos)
	% Numerical differentiation
	h = 1e-6; % since we can't divide by zero, we choice a value pretty close to 0
	J = zeros(3,3);

	% Calculate the Equations at the current position
	Eq1 = Eq{1}(pos); % equation 1
	Eq2 = Eq{2}(pos); % equation 2
	Eq3 = Eq{3}(pos); % equation 2

	%% Calculate the partial derivatives numerically
	% 1st Equation
	J(1, 1) = (Eq{1}([pos(1) + h, pos(2), pos(3)]) - Eq1) / h;
	J(1, 2) = (Eq{1}([pos(1), pos(2) + h, pos(3)]) - Eq1) / h;
	J(1, 3) = (Eq{1}([pos(1), pos(2), pos(3) + h]) - Eq1) / h;

	% 2st Equation
	J(2, 1) = (Eq{2}([pos(1) + h, pos(2), pos(3)]) - Eq2) / h;
	J(2, 2) = (Eq{2}([pos(1), pos(2) + h, pos(3)]) - Eq2) / h;
	J(2, 3) = (Eq{2}([pos(1), pos(2), pos(3) + h]) - Eq2) / h;

	% 3st Equation
	J(3, 1) = (Eq{3}([pos(1) + h, pos(2), pos(3)]) - Eq3) / h;
	J(3, 2) = (Eq{3}([pos(1), pos(2) + h, pos(3)]) - Eq3) / h;
	J(3, 3) = (Eq{3}([pos(1), pos(2), pos(3) + h]) - Eq3) / h;
end
