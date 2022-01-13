%SO3 Representation of 3D rotation
%
% This subclasss of RTBPose is an object that represents rotation in 3D.
% Internally this is a 3x3 orthonormal matrix belonging to the group SO(3).
%
% Constructor methods::
%  SO3          general constructor
%  SO3.exp      exponentiate an so(3) matrix
%  SO3.angvec   rotation about vector
%  SO3.eul      rotation defined by Euler angles
%  SO3.oa       rotation defined by o- and a-vectors
%  SO3.Rx       rotation about x-axis
%  SO3.Ry       rotation about y-axis
%  SO3.Rz       rotation about z-axis
%  SO3.rand     random orientation
%  SO3.rpy      rotation defined by roll-pitch-yaw angles
%  new          new SO3 object from instance
%
% Display and print methods::
%  plot        ^graphically display coordinate frame for pose
%  animate     ^graphically animate coordinate frame for pose
%  print       ^print the pose in single line format
%  display     ^print the pose in human readable matrix form
%
% Group operations::
%  *           ^mtimes: multiplication (group operator, transform point)
%  .*           times: multiplication (group operator) followed by normalization
%  /            ^mrdivide: multiply by inverse
%  ./           rdivide: multiply by inverse followed by normalization
%  ^            ^mpower: exponentiate (integer only)
%  .^           power: exponentiate followed by normalization
%  inv          ^inverse rotation
%  prod         ^product of elements
%
% Methods::
%  det          determinant of matrix value (is 1)
%  eig          eigenvalues of matrix value
%  interp       interpolate between rotations
%  log          logarithm of matrix value
%  norm         normalize matrix
%  simplify     ^apply symbolic simplication to all elements
%  subs         ^symbolic substitution
%  vpa          ^symbolic variable precision arithmetic
%
% Information and test methods::
%  dim         ^returns 3
%  isSE        ^returns false
%  issym       ^test if rotation matrix has symbolic elements
%  SO3.isa     test if matrix is SO(3)
%
% Conversion methods::
%  char             ^convert to human readable matrix as a string
%  SO3.convert      convert SO3 object or SO(3) matrix to SO3 object
%  double           convert to rotation matrix
%  R                convert to rotation matrix
%  SE3              convert to SE3 object with zero translation
%  T                convert to homogeneous transformation matrix with zero translation
%  toangvec         convert to rotation about vector form
%  toeul            convert to Euler angles
%  torpy            convert to roll-pitch-yaw angles
%  UnitQuaternion   convert to UnitQuaternion object
%
% Compatibility methods::
%  isrot           ^returns true
%  ishomog         ^returns false
%  trprint         ^print single line representation
%  tformplot          ^plot coordinate frame
%  tranimate       ^animate coordinate frame
%  tr2eul          convert to Euler angles
%  tr2rpy          convert to roll-pitch-yaw angles
%  tformnorm          normalize rotation matrix
%
% Operators::
%  +               ^plus: elementwise addition, result is a matrix
%  -               ^minus: elementwise subtraction, result is a matrix
%  ==              ^eq: test equality
%  ~=              ^ne: test inequality
%
% ^ inherited from RTBPose class.
%
% Properties::
%  n              normal (x) vector
%  o              orientation (y) vector
%  a              approach (z) vector
%
% See also SE2, SO2, SE3, RTBPose.

% Copyright (C) 1993-2019 Peter I. Corke
%
% This file is part of The Spatial Math Toolbox for MATLAB (SMTB).
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% https://github.com/petercorke/spatial-math

classdef SO3 < RTBPose
    
    properties (Dependent=true)
        n
        o
        a
    end
    
    methods
        
        function obj = SO3(a, varargin)
            %SO3.SO3  Construct SO3 object
            %
            % P = SO3() is the identity element, a null rotation.
            %
            % P = SO3(R) is an SO3 object formed from the rotation
            % matrix R (3x3).
            %
            % P = SO3(T) is an SO3 object formed from the rotational part
            % of the homogeneous transformation matrix T (4x4).
            %
            % P = SO3(Q) is an SO3 object that is a copy of the SO3 object Q.
            %
            % Notes::
            %  - For matrix arguments R or T the rotation submatrix is checked for validity.
            %
            % See also SE3, SO2.
            
            if nargin == 0
                % null rotation
                obj.data = eye(3,3);
            elseif isa(a, 'SO3')
                % copy
                obj.data = a.data;
            elseif SO3.isa(a)
                % from rotation matrix
                for i=1:size(a, 3)
                    x = a(:,:,i);
                    assert(SO3.isa(x, 'valid'), 'SMTB:SO3.SO3:badarg', 'matrix is not in SO(3)');
                    obj(i).data = x;
                end
            elseif SE3.isa(a)
                % from homogeneous transformation matrix, rotational part
                for i=1:size(a, 3)
                    x = a(1:3,1:3,i);
                    assert(SO3.isa(x, 'valid'), 'SMTB:SO3.SO3:badarg', 'rotation submatrix is not in SO(3)');
                    obj(i).data = x;
                end
            end
        end
        
        function printline(obj, varargin)
            trprint(double(obj), varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  GET AND SET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function RR = R(obj)
            %SO3.R  Get rotation matrix
            %
            % R = P.R() is the rotation matrix (3x3) associated with the SO3 object P.  If P
            % is a vector (1xN) then R (3x3xN) is a stack of rotation matrices, with
            % the third dimension corresponding to the index of P.
            %
            % See also SO3.T.
            if issym(obj)
                RR = zeros(3,3,length(obj), 'sym');
            else
                RR = zeros(3,3,length(obj));
            end
            for i=1:length(obj)
                RR(:,:,i) = obj(i).data(1:3,1:3);
            end
        end
        
        function TT = T(obj)
            %SO3.T  Get homogeneous transformation matrix
            %
            % T = P.T() is the homogeneous transformation matrix (4x4) associated with the
            % SO3 object P, and has zero translational component.  If P is a vector
            % (1xN) then T (4x4xN) is a stack of rotation matrices, with the third
            % dimension corresponding to the index of P.
            %
            % See also SO3.T.
            TT = zeros(4,4,length(obj));
            for i=1:length(obj)
                TT(1:3,1:3,i) = obj(i).data(1:3,1:3);
                TT(4,4,i) = 1;
            end
        end
        
        function v = get.n(obj)
            %SO3.n  Get normal vector
            %
            % P.n is the normal vector (3x1), the first column of the rotation matrix,
            % which is the x-axis unit vector.
            %
            % See also SO3.o, SO3.a.
            v = obj.data(1:3,1);
        end
        
        function v = get.o(obj)
            %SO3.o  Get orientation vector
            %
            % P.o is the orientation vector (3x1), the second column of the rotation matrix,
            % which is the y-axis unit vector..
            %
            % See also SO3.n, SO3.a.
            v = obj.data(1:3,2);
        end
        
        function v = get.a(obj)
            %SO3.a  Get approach vector
            %
            % P.a is the approach vector (3x1), the third column of the rotation matrix,
            % which is the z-axis unit vector.
            %
            % See also SO3.n, SO3.o.
            v = obj.data(1:3,3);
        end
        
        function T = norm(obj)
            %SO3.norm  Normalize rotation
            %
            % P.norm() is an SO3 object equivalent to P but with a rotation
            % matrix guaranteed to be orthogonal.
            %
            % Notes::
            %  - Overrides the classic RTB function tformnorm for an SO3 object.
            %
            % See also tformnorm.
            for k=1:length(obj)
                T(k) = SO3( tformnorm( double(obj(k))) );
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  CLASSIC RTB FUNCTION COMPATIBILITY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function T = tformnorm(obj)
            %SO3.tformnorm  Normalize rotation (compatibility)
            %
            % tformnorm(P) is an SO3 object equivalent to P but with a rotation
            % matrix guaranteed to be orthogonal.
            %
            % Notes::
            %  - Overrides the classic RTB function tformnorm for an SO3 object.
            %
            % See also tformnorm.
            for k=1:length(obj)
                T(k) = SO3( tformnorm( double(obj(k))) );
            end
        end
        
        function rpy = tr2rpy(obj, varargin)
            %SO3.tr2rpy  Convert to RPY angles (compatibility)
            %
            % tr2rpy(P, OPTIONS) is a vector (1x3) of roll-pitch-yaw angles
            % equivalent to the rotation P (SO3 object).
            %
            % Notes::
            %  - Overrides the classic RTB function tr2rpy for an SO3 object.
            %  - All the options of tr2rpy apply.
            %  - Defaults to ZYX order.
            %
            % See also tr2rpy.
            rpy = tr2rpy(obj.R, varargin{:});
        end
        
        function eul = tr2eul(obj, varargin)
            %SO3.tr2eul  Convert to Euler angles (compatibility)
            %
            % tr2eul(P, OPTIONS) is a vector (1x3) of ZYZ Euler angles
            % equivalent to the rotation P (SO3 object).
            %
            % Notes::
            %  - Overrides the classic RTB function tr2eul for an SO3 object.
            %  - All the options of tr2eul apply.
            %
            % See also tr2eul.
            eul = tr2eul(obj.R, varargin{:});
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  COMPOSITION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function out = times(obj, a)
            %SO3.times  Compose SO3 objects and normalize
            %
            % R = P1 .* P2 is an SO3 object representing the composition of the two
            % rotations described by the SO3 objects P1 and P2. This is matrix multiplication
            % of their orthonormal rotation matrices followed by normalization.
            %
            % If either, or both, of P1 or P2 are vectors, then the result is a vector.
            %  - if P1 is a vector (1xN) then R is a vector (1xN) such that R(i) = P1(i).*P2.
            %  - if P2 is a vector (1xN) then R is a vector (1xN) such that R(i) = P1.*P2(i).
            %  - if both P1 and P2 are vectors (1xN) then R is a vector (1xN) such
            %    that R(i) = P1(i).*P2(i).
            %
            % Notes::
            %  - Overloaded operator '.*'.
            %  - This is a group operator: P, Q and result all belong to the SO(3) group.
            %
            % See also RTBPose.mtimes, SO3.divide, tformnorm.
            
            % do the multiplication
            out = mtimes(obj, a);
            
            % now normalize
            if isa(out, 'SO3')
                % created an array of SE3's
                for i=1:length(out)
                    out(i).data = tformnorm(out(i).data);
                end
            end
        end
        
        
        function out = rdivide(obj, a)
            %SO3.mrdivide  Compose SO3 object with inverse and normalize
            %
            % P ./ Q is an SO3 object representing the composition of SO3 object P by the
            % inverse of SO3 object Q. This is matrix multiplication
            % of their orthonormal rotation matrices followed by normalization.
            %
            % If either, or both, of P1 or P2 are vectors, then the result is a vector.
            %  - if P1 is a vector (1xN) then R is a vector (1xN) such that R(i) = P1(i).*P2.
            %  - if P2 is a vector (1xN) then R is a vector (1xN) such that R(i) = P1.*P2(i).
            %  - if both P1 and P2 are vectors (1xN) then R is a vector (1xN) such
            %    that R(i) = P1(i).*P2(i).
            %
            % Notes::
            %  - Overloaded operator './'.
            %  - This is a group operator: P, Q and result all belong to the SO(3) group.
            %
            % See also SO3.mrdivide, SO3.times, tformnorm.
            
            % do the division
            out = mrdivide(obj, a);
            
            % now normalize
            % created an array of SO3's
            for i=1:length(out)
                out(i).data = tformnorm(out(i).data);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  OPERATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ir = inv(obj)
            %SO3.inv  Inverse
            %
            % Q = inv(P) is an SO3 object representing the inverse of the SO3 object P.
            %
            %
            % Notes::
            %  - This is a group operator: input and output in the SO(3) group.
            %  - This is simply the transpose of the underlying matrix.
            %  - P*Q will be the identity group element (zero rotation, identity matrix).
            ir = SO3(obj.R');
        end
        
        function d = det(obj)
            %SO3.inv  Determinant
            %
            % det(P) is the determinant of the SO3 object P and should always be +1.
            d = det(obj.R);
        end
        
        function R = interp(obj1, obj2, s)
            %SO3.interp Interpolate between rotations
            %
            % P1.interp(P2, s) is an SO3 object representing a slerp interpolation
            % between rotations represented by SO3 objects P1 and P2.  s varies from 0
            % (P1) to 1 (P2).  If s is a vector (1xN) then the result will be a vector
            % of SO3 objects.
            %
            % P1.interp(P2,N) as above but returns a vector (1xN) of SO3 objects
            % interpolated between P1 and P2 in N steps.
            %
            % Notes::
            %  - It is an error if any element of S is outside the interval 0 to 1.
            %
            % See also UnitQuaternion.
            
            if (length(s) == 1) && (floor(s) == s) && (s > 1)
                % is an integer, interpolate a sequence this long
                s = linspace(0, 1, s);
            else
                assert(all(s>=0 & s<=1), 'SMTB:SO3:interp:badarg', 's must be in the interval [0,1]');
            end
            q1 = quaternion(obj1.R, "rotmat", "point")
            q2 = quaternion(obj2.R, "rotmat", "point")
            qi = q1.slerp(q2, s)
            R = SO3( qi.rotmat("point") );
        end
        
        function varargout = eig(obj, varargin)
            %SO3.eig  Eigenvalues and eigenvectors
            %
            % E = eig(P) is a column vector containing the eigenvalues of the
            % underlying rotation matrix.
            %
            % [V,D] = eig(P) produces a diagonal matrix D of eigenvalues and
            % a full matrix V whose columns are the corresponding eigenvectors
            % such that A*V = V*D.
            %
            % See also eig.
            [varargout{1:nargout}] = eig(obj.data, varargin{:});
        end
        
        function S = log(obj)
            %SO3.log  Logarithm
            %
            % P.log() is the Lie algebra corresponding to the SO3 object P. It is
            % a skew-symmetric matrix (3x3).
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
            %
            % See also SO3.exp, Twist, logm, skew2vec.
            
            S = logm(obj.data);
        end
        
        
        
        %         function o = set.R(obj, data)
        %             obj.data(1:3,1:3) = data;
        %             o = obj;
        %         end
        
        %         function TT = T(obj)
        %             TT = eye(4,4);
        %             TT(1:3,1:3) = obj.data(1:3,1:3);
        %          end
        
        % TODO singularity for XYZ case,
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  conversion methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function rpy = rpy(obj, order)
            %SO3.RPY Convert to roll-pitch-yaw angles
            %
            % RPY = P.rpy(options) are the roll-pitch-yaw angles (1x3) in radians corresponding
            % to the rotational part of the SO3 object P. The 3 angles RPY=[ROLL,PITCH,YAW]
            % correspond to sequential rotations about the Z, Y and X axes
            % respectively.
            %
            % If P is a vector (1xN) then each row of RPY corresponds to an element of
            % the vector.
            %
            % Options::
            %  'xyz'   Return solution for sequential rotations about X, Y, Z axes
            %  'yxz'   Return solution for sequential rotations about Y, X, Z axes
            %
            % Notes::
            % - There is a singularity for the case where PITCH=pi/2 in which case ROLL is arbitrarily
            %   set to zero and YAW is the sum (ROLL+YAW).
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p37-38.
            %
            % See also SO3.eul, rotm2eul.
            if nargin == 1
                order = 'zyx';
            end
            rpy = rotm2eul(obj.R, order);
        end
       
        
         function rpy = rpyd(obj, order)
            %SO3.RPY Convert to roll-pitch-yaw angles
            %
            % RPY = P.rpyd(options) are the roll-pitch-yaw angles (1x3) in degrees corresponding
            % to the rotational part of the SO3 object P. The 3 angles RPY=[ROLL,PITCH,YAW]
            % correspond to sequential rotations about the Z, Y and X axes
            % respectively.
            %
            % If P is a vector (1xN) then each row of RPY corresponds to an element of
            % the vector.
            %
            % Options::
            %  'xyz'   Return solution for sequential rotations about X, Y, Z axes
            %  'yxz'   Return solution for sequential rotations about Y, X, Z axes
            %
            % Notes::
            % - There is a singularity for the case where PITCH=pi/2 in which case ROLL is arbitrarily
            %   set to zero and YAW is the sum (ROLL+YAW).
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p37-38.
            %
            % See also SO3.eul, rotm2eul.
            if nargin == 1
                order = 'zyx';
            end
            rpy = rad2deg(rotm2eul(obj.R, order));
         end
        
        %## TODO: we have inconsistent ways of handling degrees case above
        %and below
        
        function euler = eul(obj, varargin)
            %SO3.eul Convert  to Euler angles
            %
            % EUL = P.eul(OPTIONS) are the ZYZ Euler angles (1x3) corresponding to
            % the rotational part of the SO3 object P. The three angles EUL=[PHI,THETA,PSI]
            % correspond to sequential rotations about the Z, Y and Z axes
            % respectively.
            %
            % If P is a vector (1xN) then each row of EUL corresponds to an element of
            % the vector.
            %
            % Options::
            %  'deg'      Compute angles in degrees (default radians)
            %  'flip'     Choose PHI to be in quadrant 2 or 3.
            %
            % Notes::
            %  - There is a singularity when THETA=0 in which case PHI is arbitrarily
            %    set to zero and PSI is the sum (PHI+PSI).
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p36-37.
            %
            % See also SO3.rpy, rotm2eul.
            
            euler = rotm2eul(obj.R, varargin{:});
        end
       
        function aa = axang(obj)
            %SO3.axang Convert to axis-angle form
            %
            % A = P.axang(OPTIONS) is rotation expressed in terms of a
            % rotation axis axis V (1x3) and a rotation angle THETA about
            % the axis given in the form A = [V THETA]. If P is a vector
            % (1xN) then A (Nx4) is a vector of axes and angles for
            % corresponding elements of the vector, one per row.
            %
            % [V, THETA] = P.axang(OPTIONS) as above but given in two
            % components. If P is a vector (1xN) then THETA (Nx1) is a
            % vector of angles for corresponding elements of the vector and
            % V (Nx3) are the corresponding axes, one per row.
            %
            % Options::
            % 'deg'   Return angle in degrees (default radians)
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p41-42.
            %
            % See also rotm2axang.
            
            aa = rotm2axang(obj.R);
        end
        
        
        function s = to_SE3(obj)
            %SO3.SE3 Convert to SE3 object
            %
            % Q = P.to_SE3() is an SE3 object with a rotational component given by the
            % SO3 object P, and with a zero translational component.  If P is a vector
            % of SO3 objects then Q will a same length vector of SE3 objects.
            %
            % See also SE3.
            
            for i=1:length(obj)
                s(i) = SE3();
                s(i).data = [obj(i).data zeros(3,1); 0 0 0 1];
            end
        end
        
        function q = quaternion(obj)
            %SO3.quaternion Convert to quaternion object
            %
            % Q = P.quaternion() is a quaternion object of unit magnitude
            % equivalent to the rotation described by the SO3 object P. If
            % P is a vector of SO3 objects then Q will a same length vector
            % of quaternion objects.
            %
            % See also quaternion.
            R = zeros(3, 3, length(obj));
            
            % convert UQ array to 3D array of rotation matrices
            for i=1:length(obj)
                R(:,:,i) = obj(i).R;
            end
            q = quaternion(R, 'rotmat', 'point')
        end
        
        function d = dist(obj1, obj2)
            %SO3.dist Angular distance between rotations
            %
            % P.dist(Q) is the angular distance between the SO3 objects P
            % and Q.  In the range 0 to 2.
            %
            % See also quaternion.dist.
            
            q1 = quaternion(obj1.R, 'rotmat', 'point');
            q2 = quaternion(obj2.R, 'rotmat', 'point');
            d = q1.dist(q2);
        end
        
        function n = new(~, varargin)
            %SO3.new  Construct a new object of the same type
            %
            % Create a new object of the same type as the RTBPose derived instance object.
            %
            % P.new(X) creates a new object of the same type as P, by invoking the SO3 constructor on the matrix
            % X (3x3).
            %
            % P.new() as above but assumes an identity matrix.
            %
            % Notes::
            %  - Serves as a dynamic constructor.
            %  - This method is polymorphic across all RTBPose derived classes, and
            %    allows easy creation of a new object of the same class as an existing
            %    one without needing to explicitly determine its type.
            %
            % See also SE3.new, SO2.new, SE2.new.
            n = SO3(varargin{:});
        end
    end
    
    methods (Static)
        % ICANT RECALL WHY WE NEED THIS
        %         function m = ad(s)
        %             m = skew(s);
        %         end
        
        
        function obj = Rx(varargin)
            %SO3.Rx Construct SO3 from rotation about X axis
            %
            % P = SO3.Rx(THETA) is an SO3 object representing a rotation of THETA
            % radians about the x-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SO3 objects.
            %
            % P = SO3.Rx(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SO3.Ry, SO3.Rz, rotmx.
            
            theta = varargin{1};
            args = varargin(2:end);
            for i=1:length(theta)
                obj(i) = SO3( rotmx(theta(i), args{:}) );
            end
        end
        
        function obj = Ry(varargin)
            %SO3.Ry Construct SO3 from rotation about Y axis
            %
            % P = SO3.Ry(THETA) is an SO3 object representing a rotation of THETA
            % radians about the y-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SO3 objects.
            %
            % P = SO3.Ry(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SO3.Rx, SO3.Rz, rotmy.
            
            theta = varargin{1};
            args = varargin(2:end);
            for i=1:length(theta)
                obj(i) = SO3( rotmy(theta(i), args{:}) );
            end
        end
        
        function obj = Rz(varargin)
            %SO3.Rz Construct SO3 from rotation about Z axis
            %
            % P = SO3.Rz(THETA) is an SO3 object representing a rotation of THETA
            % radians about the z-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SO3 objects.
            %
            % P = SO3.Rz(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SO3.Rx, SO3.Ry, rotmz.
            
            theta = varargin{1};
            args = varargin(2:end);
            for i=1:length(theta)
                obj(i) = SO3( rotmz(theta(i), args{:}) ); %#ok<*AGROW>
            end
        end
        
        
        function R = EUL(angles, varargin)
            %SO3.eul Construct SO3 from Euler angles
            %
            % P = SO3.EUL(PHI, THETA, PSI, OPTIONS) is an SO3 object equivalent to the
            % specified Euler angles.  These correspond to rotations about the Z, Y, Z
            % axes respectively. If PHI, THETA, PSI are column vectors (Nx1) then they
            % are assumed to represent a trajectory then P is a vector (1xN) of SO3 objects.
            %
            % P = SO3.EUL(EUL, OPTIONS) as above but the Euler angles are taken from
            % consecutive columns of the passed matrix EUL = [PHI THETA PSI].  If EUL
            % is a matrix (Nx3) then they are assumed to represent a trajectory then P
            % is a vector (1xN) of SO3 objects.
            %
            % Options::
            %  'deg'      Angles are specified in degrees (default radians)
            %
            % Note::
            %   - The vectors PHI, THETA, PSI must be of the same length.
            %
            % Reference::
            %   - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p36-37.
            %
            % See also SO3.eul, SE3.RPY, eul2rotm.
            angles = flip(angles);
            if length(varargin) > 0 && strcmp(varargin(1), 'deg')
                angles = deg2rad(angles);
            end
            R = SO3( eul2rotm(angles, 'zyz') );
        end
        
        function R = RPY(varargin)
            %SO3.rpy Construct SO3 from roll-pitch-yaw angles
            %
            % P = SO3.RPY(ROLL, PITCH, YAW, OPTIONS) is an SO3 object equivalent to the
            % specified roll, pitch, yaw angles angles. These correspond to rotations
            % about the Z, Y, X axes respectively. If ROLL, PITCH, YAW are column
            % vectors (Nx1) then they are assumed to represent a trajectory then P is a
            % vector (1xN) of SO3 objects.
            %
            % P = SO3.RPY(RPY, OPTIONS) as above but the roll, pitch, yaw angles angles
            % angles are taken from consecutive columns of the passed matrix RPY =
            % [ROLL, PITCH, YAW].  If RPY is a matrix (Nx3) then they are assumed to
            % represent a trajectory and P is a vector (1xN) of SO3 objects.
            %
            % Options::
            %  'deg'   Compute angles in degrees (radians default)
            %  'xyz'   Rotations about X, Y, Z axes (for a robot gripper)
            %  'yxz'   Rotations about Y, X, Z axes (for a camera)
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p37-38
            %
            % See also SO3.rpy, SO3.EUL, SE3.RPY, eul2rotm.
            opt.unit = 'rad';
            opt.order = {'zyx', 'xyz'};
            [opt,args] = tb_optparse(opt, varargin);
            
            angles = cell2mat(args);
            angles = flip(angles);
            if strcmp(opt.unit, 'deg')
                angles = deg2rad(angles);
            end

            R = SO3( eul2rotm(angles, opt.order) );
        end
        
        function obj = Quaternion(q, frame)
            %SO3.Quaternion Construct SO3 from quaternions
            %
            % P = SO3.Quaternion(Q) is an SO3 object equivalent to the
            % quaternion array.
            %
            % P = SO3.rpy(RPY, FRAME) as above but the frame can be set
            % explicitly to either "point" (default) or "frame".
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p37-38
            %
            % See also quaternion.rotmat, SO3.EUL, SE3.RPY, TR2RPY, EUL2TR.
            if nargin < 3
                frame = "point";
            end
            obj = SO3(q.rotmat(frame));
        end
        
        function obj = OA(o, a)
            %SO3.OA Construct SO3 from orientation and approach vectors
            %
            % P = SO3.OA(O, A) is an SO3 object for the specified
            % orientation and approach vectors (3x1) formed from 3 vectors such that
            % R = [N O A] and N = O x A.
            %
            % Notes::
            %  - The rotation matrix is guaranteed to be orthonormal so long as O and A
            %   are not parallel.
            %  - The vectors O and A are parallel to the Y- and Z-axes of the coordinate
            %   frame.
            %
            % References::
            %  - Robot manipulators: mathematis, programming and control
            %    Richard Paul, MIT Press, 1981.
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p40-41.
            %
            % See also oa2rotm, SE3.OA.
            
            R = oa2rotm(o, a);
            obj = SO3(R);
        end
        
        function obj = Exp(S)
            %SO3.Exp  Construct SO3 from Lie algebra
            %
            % R = SO3.Exp(X) is the SO3 rotation corresponding to the so(3)
            % Lie algebra element SIGMA (3x3).
            %
            % R = SO3.Exp(TW) as above but the Lie algebra is represented
            % as a twist vector TW (3x1).
            %
            % Notes::
            %  - TW is the non-zero elements of X.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
            %
            % See also expm, vec2skew.
            if numel(S) == 3
                S = vec2skew(S);
            end
            obj = SO3( expm(S) );
        end
        
        function obj = Axang(theta, k)
            %SO3.angvec Construct SO3 from axis vector and angle
            %
            % R = SO3.Axang(V, THETA) is an SO3 object representing a rotation
            % of THETA about the vector V.
            %
            % R = SO3.Axang([V, THETA]) as above but construct SO3 from
            % a 4-vector.
            %
            % Notes::
            %  - If THETA == 0 then return null group element (zero rotation, identity matrix).
            %  - If THETA ~= 0 then V must have a finite length, does not have to be unit length.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p41-42.
            %
            % See also axang2rotm, SE3.Axang, eul2rotmr, rotm2axang.
            
            R = axang2rotm([k(:)' theta]);
            obj = SO3(R);
        end
        
        function T = Rand()
            %SO3.rand Construct random SO3
            %
            % SO3.Rand() is an SO3 object with a random orientation drawn from
            % a uniform distribution.
            %
            % See also randrot, quaternion.rotmat.
 
            
            T = SO3(randrot().rotmat("point"))
        end
        
        
        function R = convert(tr)
            %SO3.convert  Convert value to SO3
            %
            % Q = SO3.convert(X) is an SO3 object equivalent to X where X is either
            % an SO3 object, an SO(3) rotation matrix (3x3), an SE3 object, or an
            % SE(3) homogeneous transformation matrix (4x4).
            if isa(tr, 'SO3')
                R = SO3(tr);        % enforce it being an SO3
            elseif SO3.isa(tr)
                R = SO3(tr);
            elseif SE3.isa(tr)
                R = SO3( t2r(tr) );
            else
                error('SMTB:SO3:convert:badarg', 'expecting SO3, 3x3, SE3 or 4x4');
            end
        end
        
        function h = isa(r, ~)
            %SO3.ISA Test if a rotation matrix
            %
            % SO3.ISA(R) is true (1) if the argument is of dimension 3x3 or 3x3xN, else false (0).
            %
            % SO3.ISA(R, 'valid') as above, but also checks the validity of the rotation
            % matrix, ie. that its determinant is +1.
            %
            % Notes::
            %  - The first form is a fast, but incomplete, test for a rotation in SO(3).
            %
            % See also SE3.ISA, SE2.ISA, SO2.ISA.
            d = size(r);
            if ndims(r) >= 2
                h = all(d(1:2) == [3 3]);
                
                if h && nargin > 1 && ~isa(r, 'sym')
                    h = abs(det(r) - 1) < 10*eps;
                end
            else
                h = false;
            end
        end
        
    end
end
