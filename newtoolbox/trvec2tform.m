%TRVEC2TFORM SE(2) or SE(3) translational homogeneous transform
%
% Create a translational SE(2) matrix::
%
% T = TRANSL2(X, Y) is an SE(2) homogeneous transform (3x3) representing a
% pure translation.
%
% T = TRANSL2(P) is a homogeneous transform representing a translation or
% point P=[X,Y]. P (Mx2) represents a sequence and T (3x3xM) is a
% sequence of homogenous transforms such that T(:,:,i) corresponds to the
% i'th row of P.
%
% Extract the translational part of an SE(2) matrix::
%
% P = TRANSL2(T) is the translational part of a homogeneous transform as a
% 2-element column vector.  T (3x3xM) is a homogeneous transform
% sequence and the rows of P (Mx2) are the translational component of the
% corresponding transform in the sequence.
%
% Notes::
% - Somewhat unusually, this function performs a function and its inverse.  An
%   historical anomaly.
%
% See also SE2.t, ROTM2D, ISHOMOG2, TRPLOT2, TRANSL.

% TODO vectorize this

function T = trvec2tform(x, y, z)
    if nargin == 1
        % input is a vector of length 2 or 3
        t = x(:);
        n = length(x);

        if ~any(n == [2, 3])
            error('expecting 2 or 3 element vector')
        end
    elseif nargin == 2
        % input is a 2 values, create SE(2)
        t = [x, y]';
        n = 2;
    elseif nargin == 3
        % input is a 3 values, create SE(3)
        t = [x, y, z]';
        n = 3;
    end
    T = eye(n+1);
    T(1:n,end) = t
end
