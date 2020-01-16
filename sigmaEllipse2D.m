function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )
    %SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
    % a sigma level in the Gaussian density defined by mean and covariance.
    %
    %Input:
    %   MU          [2 x 1] Mean of the Gaussian density
    %   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
    %   LEVEL       Which sigma level curve to plot. Can take any positive value, 
    %               but common choices are 1, 2 or 3. Default = 3.
    %   NPOINTS     Number of points on the ellipse to generate. Default = 32.
    %
    %Output:
    %   XY          [2 x npoints] matrix. First row holds x-coordinates, second
    %               row holds the y-coordinates. First and last columns should 
    %               be the same point, to create a closed curve.


    %Setting default values, in case only mu and Sigma are specified.
    if nargin < 3
        level = 3;
    end
    if nargin < 4
        npoints = 32;
    end

    % Procedure:
    % - A 3 sigma level curve is given by {x} such that (x-mux)'*Q^-1*(x-mux) = 3^2
    %      or in scalar form: (x-mux) = sqrt(Q)*3
    % - replacing z= sqrtm(Q^-1)*(x-mux), such that we have now z'*z = 3^2
    %      which is now a circle with radius equal 3.
    % - Sampling in z, we have z = 3*[cos(theta); sin(theta)]', for theta=1:2*pi
    % - Back to x we get:  x = mux  + 3* sqrtm(Q)*[cos(theta); sin(theta)]'

    ang = linspace(0,2*pi,npoints);
    xy = mu + level * sqrtm(Sigma) * [cos(ang); sin(ang)];
    
    % % Alternative approach
    % [V,D] = eig(Sigma);
    % ang = linspace(0,2*pi,npoints);
    % circle = [cos(ang); sin(ang)];
    % xy = V*level*sqrt(D)*circle + mu;
end