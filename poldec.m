function [U, H, its] = poldec(A)
%POLDEC Polar decomposition
%   [U, H, its] = poldec(A) computes the polar decomposition A = U*H of the
%   square, nonsingular matrix A. ITS is the number of iterations for
%   convergence.

% Setting the parameters which can be changed if necessary
tol = 10*eps; % convergence tolerance: of order of the unit roundoff
theta = 0.6;  % switch tolerance
plots = true; % plot the errors against iterations

validateattributes(A, {'double', 'single'}, {'square'});
n = length(A);
X_prev = A;
I = eye(n);
its = 0;
terminates = false;
switches = false;
notice = false; % give a notice when switches
if plots
    fid = fopen('data.txt','w'); % store data into data.txt for plotting
end
while terminates == false
    its = its + 1;
    if ~switches % Newton method
        X_next = (X_prev + (X_prev\I)')/2;
        % compute the normwise relative distance (RD) between two
        % iterates and the equivalent forward error (FE)
        RD = norm(X_next - X_prev, Inf)/norm(X_next, Inf);
        FE = norm(I - X_next'*X_next, Inf);
        if FE <= theta
            switches = true;
            notice = true;
        end
    else % switch to Newton-Schulz method
        X_next = (3*X_prev - X_prev*(X_prev')*X_prev)/2;
        RD = norm(X_next - X_prev, Inf)/norm(X_next, Inf);
        FE = norm(I - X_next'*X_next, Inf);   
    end
    if plots
        fprintf(fid, '%.2d %.16f %.16f \n', its, RD, FE);
    end
    fprintf('%.2dth iters: Relative Distance  %.4e,  Forward Error  %.4e \n', its, RD, FE);
    if notice
        fprintf(' Now switch from the Newton iteration to the Newton-Schulz iteration.\n');
        notice = false; % only one notice is needed
    end
    if RD < tol
        terminates = true;
        if plots
            fclose(fid);
        end
    end
    X_prev = X_next;
end
U = X_prev;
H = U'*A;
H = (H' + H)/2; % Force Hermitian by taking nearest Hermitian matrix

% load data and plot the errors
if plots
    data = load('data.txt');
    semilogy(data(:,1),data(:,2),data(:,1),data(:,3));
    legend('relative distance between two iterates','equivalent forward error');
    xlabel('iterations');
    ylabel('normwise error');
end
end