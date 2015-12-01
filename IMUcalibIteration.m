%============================================================================
% Copyright (C) 2014, Heikki Hyyti
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%============================================================================

%source: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=05291740
% A Triaxial Accelerometer Calibration Method Using a Mathematical Model
% by S.P. Won, F. Golnaraghi

function [B, G] = IMUcalibIteration(S_x, S_y, S_z, B, G)
    % gravitation at Helsinki
    g0 = 9.8189; 

    % current estimate of accelerations 
    A_x = (S_x - B(1)) / G(1);
    A_y = (S_y - B(2)) / G(2);
    A_z = (S_z - B(3)) / G(3);
    
    % squares
    A_x2 = (A_x.*A_x);
    A_y2 = (A_y.*A_y);
    A_z2 = (A_z.*A_z);
    
    % error estimate
    E = A_x2 + A_y2 + A_z2 - (g0*g0);
    
    ACCEL = [A_x2 A_y2 A_z2 A_x A_y A_z];
    
    CAL = ACCEL \ E;
    
    % gain change
    G_change2 = 1 ./ (1 - CAL(1:3));
    G_change = sqrt(abs(G_change2));

    % bias change
    B_change = CAL(4:6) .* G .* G_change2 .* 0.5;
        
    % update estimates
    G = G .* G_change;
    B = B + B_change;
end