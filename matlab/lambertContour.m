function [V1, V2, info] = lambertContour(R1, R2, dt, mu, H, maxRevolutions)
    %LAMBERTCONTOURINTEGRALELLIPSE Contour solution for all conics and revolutions.
    % Five-argument calls retain their 3-by-N output. Optional maxRevolutions
    % returns 3-by-N-by-(1+2*M): zero rev, then lower-z/upper-z for each M.
    % H is a 3-by-1 reference normal (legacy) or a 1-by-N +/-1 direction vector.
    % info.status: 0 success, 1 merged roots, 5 no solution, negative failure.
    % Singular collinear geometry is unsupported: positions do not define a plane.
    % Negrete & Abdelkhalik, DOI 10.2514/1.G008499, equations 9 and 24-25.
    if nargin < 6, maxRevolutions = 0; end
    [V1, V2, info] = contourEngine(R1, R2, dt, mu, H, maxRevolutions, false);
end