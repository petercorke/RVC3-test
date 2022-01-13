%TRANSL SE(3) translational homogeneous transform
%
% Create a translational SE(3) matrix::
%
% T = TRANSL(X, Y, Z) is an SE(3) homogeneous transform (4x4) representing
% a pure translation of X, Y and Z.
%
% T = TRANSL(P) is an SE(3) homogeneous transform (4x4) representing a
% translation of P=[X,Y,Z]. P (Mx3) represents a sequence and T
% (4x4xM) is a sequence of homogeneous transforms such that T(:,:,i)
% corresponds to the i'th row of P.
%
% Extract the translational part of an SE(3) matrix::
%
% P = TRANSL(T) is the translational part of a homogeneous transform T as a
% 3-element column vector.  T (4x4xM) is a homogeneous transform
% sequence and the rows of P (Mx3) are the translational component of the
% corresponding transform in the sequence.
%
% [X,Y,Z] = TRANSL(T) is the translational part of a homogeneous transform
% T as three components.  If T (4x4xM) is a homogeneous transform sequence
% then X,Y,Z (1xM) are the translational components of the corresponding
% transform in the sequence.
%
% Notes::
% - Somewhat unusually, this function performs a function and its inverse.  An
%   historical anomaly.
%
% See also SE3.t, SE3.transl.

% TODO vectorize this


function [x, y, z] = transl(T)
    if ishomog(T)
        if ndims(T) == 3
            % transl(T)  -> P, trajectory case
            if nargout == 1
                x = squeeze(T(1:3,4,:))';
            elseif nargout == 3
                x = squeeze(T(1,4,:))';
                y = squeeze(T(2,4,:))';
                z = squeeze(T(3,4,:))';
            end
        else
            t = T(1:3, end);
            if nargout <= 1
                x = t;
            elseif nargout == 3
                x = t(1); y = t(2); z = t(3);
            else
                error('must be 1 or 3 output arguments')
            end
        end
    elseif ishomog2(T)
        if ndims(T) == 3
            % transl(T)  -> P, trajectory case
            if nargout == 1
                x = squeeze(x(1:2,3,:))';
            elseif nargout == 3
                x = squeeze(x(1,3,:))';
                y = squeeze(x(2,3,:))';
            end
        else
            t = T(1:2, end,:);
            if nargout <= 1
                x = t;
            elseif nargout == 2
                x = t(1), y = t(2)
            else
                error('must be 1 or 2 output arguments')
            end
        end
    else
        error('expecting SE(2) or SE(3) matrix')
    end
end
