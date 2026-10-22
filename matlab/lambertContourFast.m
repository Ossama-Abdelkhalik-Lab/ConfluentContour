function [V1, V2, info] = lambertContourFast(R1, R2, dt, mu, H, maxRevolutions)
    %LAMBERTCONTOURELLIPSEFAST Adaptive coarse/fine contour solver, all conics.
    % Interface and status/branch conventions match lambertContour.
    % Uses fewer initial nodes; checks residual and quadrature convergence. No
    % Newton correction, clipped roots, or silent fallback to another solver.
    if nargin < 6, maxRevolutions = 0; end
    [V1, V2, info] = contourEngine(R1, R2, dt, mu, H, maxRevolutions, true);
end