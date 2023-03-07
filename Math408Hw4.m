
%% Gradient Method with backtracking stepsize rule
% Define the problem
A = hilb(5);
x0 = [1; 2; 3; 4; 5];
tol = 1e-4;

% Set the parameters for the backtracking line search
alpha = 0.5;
beta = 0.5;
s = 1;

% Initialize the algorithm
x = x0;
grad = A * x;
iter_count = 0;

% Perform gradient descent with backtracking line search
while norm(grad) > tol
    % Determine the step size using backtracking line search
    step_size = s;
    while A * ((x - step_size * grad)' * grad) > (A * x)' * grad - alpha * step_size * norm(grad)^2
        step_size = beta * step_size;
    end
    
    % Update the solution
    x = x - step_size * grad;
    grad = A * x;
    iter_count = iter_count + 1;
end

% Display the results
disp('Solution using gradient descent with backtracking line search:');
disp(x);
fprintf('Number of iterations: %d\n', iter_count);

%% Gradient method with exact line search
% Define the problem
A = hilb(5);
x0 = [1; 2; 3; 4; 5];
tol = 1e-4;

% Initialize the algorithm
x = x0;
grad = A * x;
iter_count = 0;

% Perform gradient descent with exact line search
while norm(grad) > tol
    % Determine the step size using exact line search
    step_size = norm(grad)^2 / (grad' * A * grad);
    
    % Update the solution
    x = x - step_size * grad;
    grad = A * x;
    iter_count = iter_count + 1;
end

% Display the results
disp('Solution using gradient descent with exact line search:');
disp(x);
fprintf('Number of iterations: %d\n', iter_count);


%% Diagonally scaled gradient method with exact line search
% Define the problem
A = hilb(5);
x0 = [1; 2; 3; 4; 5];
tol = 1e-4;

% Calculate the diagonal scaling matrix
D = diag(1./diag(A));

% Initialize the algorithm
x = x0;
grad = A * x;
iter_count = 0;

% Perform diagonally scaled gradient descent with exact line search
while norm(grad) > tol
    % Determine the step size using exact line search
    grad_scaled = D * grad;
    step_size = norm(grad_scaled)^2 / (grad_scaled' * A * grad_scaled);
    
    % Update the solution
    x = x - step_size * grad_scaled;
    grad = A * x;
    iter_count = iter_count + 1;
end

% Display the results
disp('Solution using diagonally scaled gradient descent with exact line search:');
disp(x);
fprintf('Number of iterations: %d\n', iter_count);

%% Diagonally scaled gradient method with backtracking line search
% Define the problem
A = hilb(5);
x0 = [1; 2; 3; 4; 5];
tol = 1e-4;
alpha = 0.1;
beta = 0.5;
s = 1;

% Calculate the diagonal scaling matrix
D = diag(1./diag(A));

% Initialize the algorithm
x = x0;
grad = A * x;
iter_count = 0;

% Perform diagonally scaled gradient descent with backtracking line search
while norm(grad) > tol
    % Determine the step size using backtracking line search
    grad_scaled = D * grad;
    step_size = s;
    while x - step_size * grad_scaled <= 0
        step_size = beta * step_size;
    end
    while norm(A * (x - step_size * grad_scaled) - A * x) > alpha * step_size * norm(grad_scaled)^2
        step_size = beta * step_size;
    end
    
    % Update the solution
    x = x - step_size * grad_scaled;
    grad = A * x;
    iter_count = iter_count + 1;
end

% Display the results
disp('Solution using diagonally scaled gradient descent with backtracking line search:');
disp(x);
fprintf('Number of iterations: %d\n', iter_count);
