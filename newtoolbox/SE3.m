%SE3 Representation of 3D rigid-body motion
%
% This subclasss of RTBPose is an object that represents rigid-body motion in 2D.
% Internally this is a 3x3 homogeneous transformation matrix (4x4) belonging to
% the group SE(3).
%
% Constructor methods::
%  SE3              general constructor
%  SE3.angvec       rotation about vector
%  SE3.eul          rotation defined by Euler angles
%  SE3.exp          exponentiate an se(3) matrix
%  SE3.oa           rotation defined by o- and a-vectors
%  SE3.Rx           rotation about x-axis
%  SE3.Ry           rotation about y-axis
%  SE3.Rz           rotation about z-axis
%  SE3.rand         random transformation
%  SE3.rpy          rotation defined by roll-pitch-yaw angles
%  new              new SE3 object
%
% Display and print methods::
%  animate          ^graphically animate coordinate frame for pose
%  display          ^print the pose in human readable matrix form
%  plot             ^graphically display coordinate frame for pose
%  print            ^print the pose in single line format
%
% Group operations::
%  *               ^mtimes: multiplication (group operator, transform point)
%  .*              ^^times: multiplication (group operator) followed by normalization
%  /               ^mrdivide: multiply by inverse
%  ./              ^^rdivide: multiply by inverse followed by normalization
%  ^               ^mpower: xponentiate (integer only)
%  .^              ^power: exponentiate followed by normalization
%  inv             inverse
%  prod            ^product of elements
%
% Methods::
%  det              determinant of matrix component
%  eig              eigenvalues of matrix component
%  log              logarithm of rotation matrixr>=0 && r<=1ub
%  simplify         ^apply symbolic simplication to all elements
%  Ad               adjoint matrix (6x6)
%  increment        update pose based on incremental motion
%  interp           interpolate poses
%  velxform         compute velocity transformation
%  interp           interpolate between poses
%  ctraj            Cartesian motion
%  norm             normalize the rotation submatrix
%
% Information and test methods::
%  dim*             returns 4
%  isSE*            returns true
%  issym*           test if rotation matrix has symbolic elements
%  isidentity       test for null motion
%  SE3.isa          check if matrix is SE(3)
%
% Conversion methods::
%  char             convert to human readable matrix as a string
%  SE3.convert      convert SE3 object or SE(3) matrix to SE3 object
%  double           convert to SE(3) matrix
%  R                convert rotation part to SO(3) matrix
%  SO3              convert rotation part to SO3 object
%  T                convert to SE(3) matrix
%  t                translation column vector
%  toangvec         convert to rotation about vector form
%  todelta          convert to differential motion vector
%  toeul            convert to Euler angles
%  torpy            convert to roll-pitch-yaw angles
%  tv               translation column vector for vector of SE3
%  UnitQuaternion   convert to UnitQuaternion object
%
% Compatibility methods::
%  homtrans         apply to vector
%  isrot            ^returns false
%  ishomog          ^returns true
%  t2r              ^convert to rotation matrix
%  tr2rt            ^convert to rotation matrix and translation vector
%  tr2eul           ^^convert to Euler angles
%  tr2rpy           ^^convert to roll-pitch-yaw angles
%  tranimate        ^animate coordinate frame
%  transl           translation as a row vector
%  tformnorm           ^^normalize the rotation matrix
%  tformplot           ^plot coordinate frame
%  trprint          ^print single line representation
%
% Other operators::
%  +                ^plus: elementwise addition, result is a matrix
%  -                ^minus: elementwise subtraction, result is a matrix
%  ==               ^eq: test equality
%  ~=               ^ne: test inequality
%
% - ^ inherited from RTBPose
% - ^^ inherited from SO3
%
% Properties::
%  n              get.n: normal (x) vector
%  o              get.o: orientation (y) vector
%  a              get.a: approach (z) vector
%  t              get.t: translation vector
%
% For single SE3 objects only, for a vector of SE3 objects use the
% equivalent methods
% t       translation as a 3x1 vector (read/write)
% R       rotation as a 3x3 matrix (read)
%
% Notes::
%  - The properies R, t are implemented as MATLAB dependent properties.
%    When applied to a vector of SE3 object the result is a comma-separated
%    list which can be converted to a matrix by enclosing it in square
%    brackets, eg [T.t] or more conveniently using the method T.transl
%
% See also SO3, SE2, RTBPose.

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

%TODO
% interp
% animate
% chain
% animate
% jacobian
% log
% norm
% print
%
% all vectorised!
%
%
% Superclass, provides char,display, get/set R,t,T

classdef SE3 < SO3
    
    properties (Dependent = true)
        t
        %         T
        %R
    end
    
    methods
        
        function obj = SE3(varargin)
            %SE3.SE3 Create an SE(3) object
            %
            % Constructs an SE(3) pose object that contains a 4x4 homogeneous transformation matrix.
            %
            % T = SE3() is the identity element, a null motion.
            %
            % T = SE3(X, Y, Z) is an object representing pure translation defined by X,
            % Y and Z.
            %
            % T = SE3(XYZ) is an object representing pure translation defined by XYZ
            % (3x1).  If XYZ (Nx3) returns an array of SE3 objects, corresponding to
            % the rows of XYZ.
            %
            % T = SE3(T) is an object representing translation and rotation defined by
            % the homogeneous transformation matrix T (3x3).  If T (3x3xN) returns an array of SE3 objects,
            % corresponding to the third index of T.
            %
            % T = SE3(R, XYZ) is an object representing rotation defined by the
            % orthonormal rotation matrix R (3x3) and position given by XYZ (3x1).
            %
            % T = SE3(T) is a copy of the SE3 object T. If T (Nx1) returns an array of SE3 objects,
            % corresponding to the index of T.
            %
            % Options::
            % 'deg'         Angle is specified in degrees
            %
            % Notes::
            %  - Arguments can be symbolic.
            %  - R and T are checked to be valid SO(2) or SE(2) matrices.
            
            obj.data = eye(4,4);
            args = varargin;
            
            % if any of the arguments is symbolic the result will be symbolic
            if any( cellfun(@(x) isa(x, 'sym'), args) )
                obj.data = sym(obj.data);
            end
            
            switch length(args)
                case 0
                    % ()
                    obj.data = eye(4,4);
                    return;
                    
                case 1
                    a = args{1};
                    if isvec(a, 3)
                        % (t)
                        obj.t = a(:);
                    elseif SE3.isa(a)
                        % (SE3)
                        for i=1:size(a, 3)
                            obj(i).data = a(:,:,i); %#ok<*AGROW>
                        end
                    elseif isa(a, 'SE3')
                        % (T)
                        for i=1:length(a)
                            obj(i).data = a(i).data;
                        end
                    elseif isa(a, 'SO3')
                        % (R)
                        for i=1:length(a)
                            obj(i).data(1:3,1:3) = a.R;
                        end
                    elseif size(a,2) == 3
                        % SE3( xyz )
                        for i=1:length(a)
                            %obj(i).data = SE3(a(i,:));
                            obj(i).data(1:3,4) = a(i,:)';
                        end
                    else
                        error('SMTB:SE3:badarg', 'unknown arguments');
                    end
                    
                case 2
                    a = args{1}; b = args{2};
                    
                    if (isrot(a) || SO3.isa(a)) && isvec(b,3)
                        % (R, t)
                        if isrot(a)
                            obj.data(1:3,1:3) = a;
                        else
                            obj.R(1:3,1:3) = a.R;
                        end
                        obj.data(1:3,4) = b(:);
                    else
                        error('SMTB:SE3:badarg', 'unknown arguments');
                    end
                    
                case 3
                    a = args{1}; b = args{2}; c = args{3};
                    
                    obj.data(1:3,4) = [a; b; c];
                    
                otherwise
                    error('SMTB:SE3:badarg', 'too many arguments');
                    
            end
            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  GET AND SET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function o = get.t(obj)
            o = obj.data(1:3,4);
        end
        
        function o = set.t(obj, t) %#ok<MCSNOV>
            %SE3.t Get translation vector
            %
            % T = P.t is the translational part of SE3 object as a 3-element column
            % vector.
            %
            % Notes::
            %  - If applied to  a vector will return a comma-separated list, use
            %    .tv() instead.
            %
            % See also SE3.tv, transl.
            
            % TODO CHECKING
            obj.data(1:3,4) = t(:);
            o = obj;
        end
        
        
        function [tx,ty,tz] = transl(obj)
            %SE3.transl Get translation vector
            %
            % T = P.transl() is the translational part of SE3 object as a 3-element row
            % vector.  If P is a vector (1xN) then
            %  the rows of T (Mx3) are the translational component of the
            % corresponding pose in the sequence.
            %
            % [X,Y,Z] = P.transl() as above but the translational part is returned as
            % three components.  If P is a vector (1xN) then X,Y,Z (1xN) are the
            % translational components of the corresponding pose in the sequence.
            %
            % Notes::
            %  - The .t method only works for a single pose object, on  a vector it
            %    returns a comma-separated list.
            %
            % See also SE3.t, transl.
            
            if nargout == 1 || nargout == 0
                tx = [obj.t]';
            else
                t = obj.t; %#ok<*PROP>
                tx = t(1);
                ty = t(2);
                tz = t(3);
            end
        end
        
        
        function TT = T(obj)
            %SO2.T  Get homogeneous transformation matrix
            %
            % T = P.T() is the homogeneous transformation matrix (3x3) associated with the
            % SO2 object P, and has zero translational component.  If P is a vector
            % (1xN) then T (3x3xN) is a stack of rotation matrices, with the third
            % dimension corresponding to the index of P.
            %
            % See also SO2.T.
            TT = double(obj);
        end
        
        
        % set.R
        
        % handle assignment SE3 = SE3, SE3=4x4, can'toverload equals
        
        function t = tv(obj)
            %SE.tv Return translation for a vector of SE3 objects
            %
            % P.tv is a column vector (3x1) representing the translational part of the
            % SE3 pose object P.  If P is a vector of SE3 objects (Nx1) then the result
            % is a matrix (3xN) with columns corresponding to the elements of P.
            %
            % See also SE3.t.
            
            %t = zeros(3,length(obj)); FAILS FOR SYMBOLIC CASE
            for i=1:length(obj)
                t(:,i) = obj(i).data(1:3,4);
            end
        end
        
        function v = isidentity(obj)
            %SE3.identity  Test if identity element
            %
            % P.isidentity() is true if the SE3 object P corresponds to null motion,
            % that is, its homogeneous transformation matrix is identity.
            
            v = all(all(obj.T == eye(4,4)));
        end
        
        function T = increment(obj, v)
            %SE3.increment  Apply incremental motion to an SE3 pose
            %
            % P1 = P.increment(D) is an SE3 pose object formed by compounding the
            % SE3 pose with the incremental motion described by D (6x1).
            %
            % The vector D=(dx, dy, dz, dRx, dRy, dRz) represents infinitessimal translation
            % and rotation, and is an approximation to the instantaneous spatial velocity
            % multiplied by time step.
            %
            % See also SE3.todelta, SE3.delta, delta2tform, tform2delta.
            T = obj .* SE3(delta2tform(v));
        end
        
        function J = velxform(obj, varargin)
            %SE3.velxform  Velocity transformation
            %
            % Transform velocity between frames.  A is the world frame, B is the body
            % frame and C is another frame attached to the body.  PAB is the pose of
            % the body frame with respect to the world frame, PCB is the pose of the
            % body frame with respect to frame C.
            %
            % J = PAB.velxform() is a 6x6 Jacobian matrix that maps velocity from frame
            % B to frame A.
            %
            % J = PCB.velxform('samebody') is a 6x6 Jacobian matrix that maps velocity
            % from frame C to frame B.  This is also the adjoint of PCB.
            
            opt.samebody = false;
            
            opt = tb_optparse(opt, varargin);
            
            R = obj.R;
            
            if opt.samebody
                J = [
                    R           vec2skew(obj.t)*R
                    zeros(3,3)  R
                    ];
            else
                J = [
                    R            zeros(3,3)
                    zeros(3,3)   R
                    ];
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  STANDARD SE(3) MATRIX OPERATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function it = inv(obj)
            %SE3.inv  Inverse of SE3 object
            %
            % Q = inv(P) is the inverse of the SE3 object P.
            %
            % Notes::
            %  - This is formed explicitly, no matrix inverse required.
            %  - This is a group operator: input and output in the SE(3)) group.
            %  - P*Q will be the identity group element (zero motion, identity matrix).
            
            
            it = SE3( obj.R', -obj.R'*obj.t);
        end
        
        function Ti = interp(obj1, varargin)
            %SE3.interp Interpolate SE3 poses
            %
            % P1.interp(P2, s) is an SE3 object representing an interpolation
            % between poses represented by SE3 objects P1 and P2.  s varies from 0
            % (P1) to 1 (P2).  If s is a vector (1xN) then the result will be a vector
            % of SO3 objects.
            %
            % P1.interp(P2, N) as above but returns a vector (1xN) of SE3 objects
            % interpolated between P1 and P2 in N steps.
            %
            % Notes::
            %  - The rotational interpolation (slerp) can be interpretted
            %   as interpolation along a great circle arc on a sphere.
            %  - It is an error if any element of S is outside the interval 0 to 1.
            %
            % See also TRINTERP, CTRAJ, UnitQuaternion.
            
            % TODO replace with transformtraj
            
            if isa(varargin{1}, 'SE3')
                % interp(SE3, SE3, s)  interpolate between given values
                obj2 = varargin{1};
                varargin = varargin(2:end);
                try
                    Ti = SE3( trinterp(obj1.T, obj2.T, varargin{:}) );
                catch me
                    switch me.identifier
                        case 'SMTB:trinterp:badarg'
                            throw( MException('SMTB:SE3:interp:badarg', 'value of S outside interval [0,1]') );
                        otherwise
                            rethrow(me);
                    end
                end
            else
                % interp(SE3, s)  interpolate between null and given value
                try
                    T = trinterp(obj1.T, varargin{:});
                catch me
                    switch me.identifier
                        case 'SMTB:trinterp:badarg'
                            throw( MException('SMTB:SE3:interp:badarg', 'value of S outside interval [0,1]') );
                        otherwise
                            rethrow(me);
                    end
                end
                Ti = SE3(T);
            end
        end
        
        function traj = ctraj(T0, T1, t)
            %SE3.ctraj Cartesian trajectory between two poses
            %
            % TC = T0.ctraj(T1, N) is a Cartesian trajectory defined by a vector of SE3
            % objects (1xN) from pose T0 to T1, both described by SE3 objects.  There
            % are N points on the trajectory that follow a trapezoidal velocity profile
            % along the trajectory.
            %
            % TC = CTRAJ(T0, T1, S) as above but the elements of S (Nx1) specify the
            % fractional distance  along the path, and these values are in the range [0 1].
            % The i'th point corresponds to a distance S(i) along the path.
            %
            % Notes::
            %  - In the second case S could be generated by a scalar trajectory generator
            %    such as TPOLY or LSPB (default).
            %  - Orientation interpolation is performed using quaternion interpolation.
            %
            % Reference::
            % Robotics, Vision & Control, Sec 3.1.5,
            % Peter Corke, Springer 2011
            %
            % See also LSPB, MSTRAJ, TRINTERP, CTRAJ, UnitQuaternion.interp.
            
            assert(isa(T1, 'SE3'), 'second transform must be SE3');
            
            % distance along path is a smooth function of time
            if isscalar(t)
                s = lspb(0, 1, t);
            else
                s = t(:);
            end
            
            for i=1:length(s)
                traj(i) = T0.interp(T1, s(i));
            end
        end
        
        function m = log(obj)
            %SE3.log  Lie algebra
            %
            % P.log() is the Lie algebra corresponding to the SE3 object P. It is
            % an augmented skew-symmetric matrix (4x4).
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
            %
            % See also SE3.logs, SE3.Twist, logm.
            m = logm(obj.T);
        end
        
        function m = logs(obj)
            %SE3.log  Lie algebra in vector form
            %
            % P.logs() is the Lie algebra expressed as a vector (1x6)
            % corresponding to the SE3 object P.  The vector comprises the
            % translational elements followed by the unique elements of the
            % skew-symmetric upper-left 3x3 submatrix.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
            %
            % See also SE3.log, SE3.Twist, logm, skewa2vec.
            
            m = logm(obj.T);
            m = skewa2vec(m);
        end
        
        function tw = Twist(obj)
            %SE3.Twist  Convert to Twist object
            %
            % TW = P.Twist() is the equivalent Twist object.  The elements of the twist are the unique
            % elements of the Lie algebra of the SE3 object P.
            %
            % See also SE3.logs, Twist.
            tw = Twist( obj.log );
        end
        
        
        function m = adjoint(obj)
            %SE3.Ad  Adjoint matrix
            %
            % A = P.Ad() is the adjoint matrix (6x6) corresponding to the pose P.
            %
            % See also Twist.ad.
            m = [
                obj.R       vec2skew(obj.t)*obj.R
                zeros(3,3)  obj.R
                ];
        end
        
        %% ad function
        %         function m = Ad(obj)
        %             m = obj.R;
        %         end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Rotation conversion wrappers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        function out = eul(obj, varargin)
            %SE3.eul Convert  to Euler angles
            %
            % EUL = P.eul(OPTIONS) are the ZYZ Euler angles (1x3) corresponding to
            % the rotational part of the SE3 object P. The 3 angles EUL=[PHI,THETA,PSI]
            % correspond to sequential rotations about the Z, Y and Z axes
            % respectively.
            %
            % If P is a vector (1xN) then each row of EUL corresponds to an element of
            % the vector.
            %
            % Options::
            %  'deg'      Compute angles in degrees (radians default)
            %  'flip'     Choose first Euler angle to be in quadrant 2 or 3.
            %
            % Notes::
            %  - There is a singularity for the case where THETA=0 in which case PHI is arbitrarily
            %    set to zero and PSI is the sum (PHI+PSI).
            %
            % See also SO3.eul, SE3.rpy, tform2eul.
            out = obj.to_SO3.eul(varargin{:});
        end
        
        function out = rpy(obj, varargin)
            %SE3.RPY Convert to roll-pitch-yaw angles
            %
            % RPY = P.rpy(options) are the roll-pitch-yaw angles (1x3) corresponding
            % to the rotational part of the SE3 object P. The 3 angles RPY=[R,P,Y]
            % correspond to sequential rotations about the Z, Y and X axes
            % respectively.
            %
            % If P is a vector (1xN) then each row of RPY corresponds to an element of
            % the vector.
            %
            % Options::
            %  'deg'   Compute angles in degrees (radians default)
            %  'xyz'   Return solution for sequential rotations about X, Y, Z axes
            %  'yxz'   Return solution for sequential rotations about Y, X, Z axes
            %
            % Notes::
            %  - There is a singularity for the case where P=pi/2 in which case R is arbitrarily
            %    set to zero and Y is the sum (R+Y).
            %
            % See also SO3.rpy, SE3.eul, tform2eul.
            out = obj.to_SO3.rpy(varargin{:});
        end
        
        function varargout = axang(obj, varargin)
            %SE3.axang Convert to axis-angle form
            %
            % A = P.axang(OPTIONS) is rotational part of P expressed in terms of a
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
            % See also tform2axang.
            [varargout{1:nargout}] = obj.to_SO3.toangvec(varargin{:});

        end
        
        function d = delta(P1, P2)
            %SE3.delta Convert SE3 object to differential motion vector
            %
            % D = P0.delta(P1) is the differential motion (6x1) corresponding to
            % infinitessimal motion (in the P0 frame) from SE3 pose P0 to P1.
            %
            % The vector D=(dx, dy, dz, dRx, dRy, dRz) represents infinitessimal translation
            % and rotation, and is an approximation to the instantaneous spatial velocity
            % multiplied by time step.
            %
            % D = P.delta() as above but the motion is from the world frame to the SE3
            % pose P.
            %
            % Notes::
            %  - D is only an approximation to the motion, and assumes
            %    that P0 ~ P1 or P ~ eye(4,4).
            %  - can be considered as an approximation to the effect of spatial velocity over a
            %    a time interval, average spatial velocity multiplied by time.
            %
            % See also SE3.increment, tform2delta, delta2tform.
            
            if nargin == 1
                d = tform2delta(P1.T);
            elseif nargin == 2
                d = tform2delta(P1.T, P2.T);
            end
        end
        
        % conversion methods
        
        function s = to_SO3(obj)
            %SE3.to_SO3 Convert rotational component to SO3 object
            %
            % P.to_SO3 is an SO3 object representing the rotational component of the SE3
            % pose P.  If P is a vector (Nx1) then the result is a vector (Nx1).
            for i=1:length(obj)
                s(i) = SO3();
                s(i).data= obj(i).R;
            end
        end
        
        function n = new(~, varargin)
            %SE3.new  Construct a new object of the same type
            %
            % P2 = P.new(X) creates a new object of the same type as P, by invoking the SE3 constructor on the matrix
            % X (4x4).
            %
            % P2 = P.new() as above but defines a null motion.
            %
            % Notes::
            %  - Serves as a dynamic constructor.
            %  - This method is polymorphic across all RTBPose derived classes, and
            %    allows easy creation of a new object of the same class as an existing
            %    one without needing to explicitly determine its type.
            %
            % See also SO3.new, SO2.new, SE2.new.
            
            n = SE3(varargin{:});
        end
        
        function T = normalize(obj)
            %SE3.normalize  Normalize rotation submatrix (compatibility)
            %
            % P.normalize() is an SE3 pose equivalent to P but the rotation
            % matrix is normalized (guaranteed to be orthogonal).
            %
            % See also tformnorm.
            for k=1:length(obj)
                T(k) = SE3( tformnorm( double(obj(k))) );
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  VANILLA RTB FUNCTION COMPATIBILITY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function T = tformnorm(obj)
            %SE3.tformnorm  Normalize rotation submatrix (compatibility)
            %
            % T = tformnorm(P) is an SE3 object equivalent to P but
            % normalized (rotation matrix guaranteed to be orthogonal).
            %
            % Notes::
            %  - Overrides the classic RTB function tformnorm for an SE3 object.
            %
            % See also tformnorm.
            for k=1:length(obj)
                T(k) = SE3( tformnorm( double(obj(k))) );
            end
        end
        
        function Vt = homtrans(P, V)
            %SE3.homtrans  Apply transformation to points (compatibility)
            %
            % homtrans(P, V) applies SE3 pose object P to the points stored columnwise in
            % V (3xN) and returns transformed points (3xN).
            %
            % Notes::
            %  - P is an SE3 object defining the pose of {A} with respect to {B}.
            %  - The points are defined with respect to frame {A} and are transformed to be
            %    with respect to frame {B}.
            %  - Equivalent to P*V using overloaded SE3 operators.
            %
            % See also RTBPose.mtimes, HOMTRANS.
            Vt = P*V;
        end
        
        
    end
    
    methods (Static)
        
        
        % Static factory methods for constructors from standard representations
        
        function obj = Rx(varargin)
            %SE3.Rx Construct SE3 from rotation about X axis
            %
            % P = SE3.Rx(THETA) is an SE3 object representing a rotation of THETA
            % radians about the x-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SE3 objects.
            %
            % P = SE3.Rx(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SE3.Ry, SE3.Rz, rotmx.
            obj = SE3( SO3.Rx(varargin{:}) );
        end
        
        function obj = Ry(varargin)
            %SE3.Ry Construct SE3 from rotation about Y axis
            %
            % P = SE3.Ry(THETA) is an SE3 object representing a rotation of THETA
            % radians about the y-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SE3 objects.
            %
            % P = SE3.Ry(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SE3.Ry, SE3.Rz, rotmy.
            obj = SE3( SO3.Ry(varargin{:}) );
        end
        
        function obj = Rz(varargin)
            %SE3.Rz Construct SE3 from rotation about Z axis
            %
            % P = SE3.Rz(THETA) is an SE3 object representing a rotation of THETA
            % radians about the z-axis.  If the THETA is a vector (1xN) then P will be
            % a vector (1xN) of corresponding SE3 objects.
            %
            % P = SE3.Rz(THETA, 'deg') as above but THETA is in degrees.
            %
            % See also SE3.Ry, SE3.Rz, rotmz.
            obj = SE3( SO3.Rz(varargin{:}) );
        end
        
        function obj = OA(varargin)
            %SE3.OA Construct SE3 from orientation and approach vectors
            %
            % P = SE3.OA(O, A) is an SE3 object for the specified
            % orientation and approach vectors (3x1) formed from 3 vectors such that R
            % = [N O A] and N = O x A, with zero translation.
            %
            % Notes::
            %  - The rotation submatrix is guaranteed to be orthonormal so long as O and A
            %    are not parallel.
            %  - The vectors O and A are parallel to the Y- and Z-axes of the coordinate
            %    frame.
            %
            % References::
            %  - Robot manipulators: mathematics, programming and control
            %    Richard Paul, MIT Press, 1981.
            %
            % See also oa2tform, SO3.oa.

            T = oa2tform(o, a);
            obj = SE3(R);
        end
        
        function obj = Axang(varargin)
            %SE3.Axang Construct SE3 from axis vector and angle
            %
            % T = SO3.Axang(V, THETA) is an SE3 object representing a rotation
            % of THETA about the vector V.
            %
            % R = SO3.Axang([V, THETA]) as above but construct SE3 from
            % a 4-vector.
            %
            % Notes::
            %  - If THETA == 0 then return null group element (zero rotation, identity matrix).
            %  - If THETA ~= 0 then V must have a finite length, does not have to be unit length.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p41-42.
            %
            % See also axang2tform, SO3.Axang, eul2tform, tform2axang.
            obj = SE3( SO3.angvec(varargin{:}) );
        end
        
        function obj = RPY(varargin)
            %SE3.rpy Construct SE3 from roll-pitch-yaw angles
            %
            % P = SE3.RPY(ROLL, PITCH, YAW, OPTIONS) is an SE3 object equivalent to the
            % specified roll, pitch, yaw angles angles with zero translation. These correspond to rotations
            % about the Z, Y, X axes respectively. If ROLL, PITCH, YAW are column
            % vectors (Nx1) then they are assumed to represent a trajectory then P is a
            % vector (1xN) of SE3 objects.
            %
            % P = SE3.RPY(RPY, OPTIONS) as above but the roll, pitch, yaw angles angles
            % angles are taken from consecutive columns of the passed matrix RPY =
            % [ROLL, PITCH, YAW].  If RPY is a matrix (Nx3) then they are assumed to
            % represent a trajectory and P is a vector (1xN) of SE3 objects.
            %
            % Options::
            %  'deg'   Compute angles in degrees (radians default)
            %  'xyz'   Rotations about X, Y, Z axes (for a robot gripper)
            %  'yxz'   Rotations about Y, X, Z axes (for a camera)
            %
            % Reference::
            % - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p37-38.
            %
            % See also SE3.rpy, SE3.EUL, eul2tform.
            obj = SE3( SO3.RPY(varargin{:}) );
        end
        
        function obj = EUL(varargin)
            %SE3.EUL Construct SE3 from Euler angles
            %
            % P = SO3.EUL(PHI, THETA, PSI, OPTIONS) is an SE3 object equivalent to the
            % specified Euler angles.  These correspond to rotations about the Z, Y, Z
            % axes respectively. If PHI, THETA, PSI are column vectors (Nx1) then they
            % are assumed to represent a trajectory then P is a vector (1xN) of SE3 objects.
            %
            % P = SO3.EUL(EUL, OPTIONS) as above but the Euler angles are taken from
            % consecutive columns of the passed matrix EUL = [PHI THETA PSI].  If EUL
            % is a matrix (Nx3) then they are assumed to represent a trajectory then P
            % is a vector (1xN) of SE3 objects.
            %
            % Options::
            %  'deg'      Angles are specified in degrees (default radians)
            %
            % Note::
            %  - Translation is zero.
            %  - The vectors PHI, THETA, PSI must be of the same length.
            %
            % Reference::
            % - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p36-37.
            %
            % See also SO3.eul, SE3.RPY, eul2tform.
            obj = SE3( SO3.EUL(varargin{:}) );
        end
        
        function T = Rand()
            %SE3.Rand Construct random SE3
            %
            % SE3.Rand() is an SE3 object with a uniform random translation and a
            % uniform random RPY/ZYX orientation.  Translation components are in the interval -1 to
            % 1.
            %
            % See also randrot, quaternion.rotmat.
            T = SE3( 2*rand(3,1)-1 ) * SE3(SO3.Rand);
        end
        
        function obj = Delta(d)
            %SE3.Delta Construct SE3 object from differential motion vector
            %
            % T = SE3.Delta(D) is an SE3 pose object representing differential
            % motion D (6x1).
            %
            % The vector D=(dx, dy, dz, dRx, dRy, dRz) represents infinitessimal translation
            % and rotation, and is an approximation to the instantaneous spatial velocity
            % multiplied by time step.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p67.
            %
            % See also SE3.delta, SE3.increment, tform2delta.
            
            assert(isvec(d,6), 'SMTB:SE3:badarg', 'delta is a 6-vector');
            obj = SE3( delta2tform(d));
        end
        
        % Static factory methods for constructors from exotic representations
        function obj = Exp(S)
            %SE3.exp Construct SE3 from Lie algebra
            %
            % SE3.exp(SIGMA) is the SE3 rigid-body motion corresponding to the se(3)
            % Lie algebra element SIGMA (4x4).
            %
            % SE3.exp(TW) as above but the Lie algebra is represented
            % as a twist vector TW (6x1).
            %
            % SE3.exp(SIGMA, THETA) as above, but the motion is given by SIGMA*THETA
            % where SIGMA is an se(3) element (4x4) whose rotation part has a unit norm.
            %
            % Notes::
            %  - TW is the non-zero elements of X.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
            %
            % See also trexp, vec2skewa, Twist.
            
            % TODO deal with all cases
            if numel(S) == 6
                S = vec2skewa(S);
            end
            obj = SE3( expm(S) );
        end
        
        %         function m = ad(s)
        %             %m = [skew(s(1:3)) skew(s(4:6));  skew(s(4:6)) zeros(3,3)];
        %             m = [skew(s(4:6)) skew(s(1:3)) ;  zeros(3,3) skew(s(4:6)) ];
        %         end
        
        
        function T = convert(tr)
            %SE3.check  Convert to SE3
            %
            % Q = SE3.convert(X) is an SE3 object equivalent to X where X is either
            % an SE3 object, or an SE(3) homogeneous transformation matrix (4x4).
            if isa(tr, 'SE3')
                T = tr;
            elseif SE3.isa(tr)
                T = SE3(tr);
            else
                error('expecting an SE3 or 4x4 matrix');
            end
        end
        
        function h = isa(tr, ~)
            %SE3.ISA Test if matrix is SE(3)
            %
            % SE3.ISA(T) is true (1) if the argument T is of dimension 4x4 or 4x4xN, else
            % false (0).
            %
            % SE3.ISA(T, 'valid') as above, but also checks the validity of the rotation
            % sub-matrix.
            %
            % Notes::
            %  - Is a class method.
            %  - The first form is a fast, but incomplete, test for a transform in SE(3).
            %
            % See also SO3.isa, SE2.isa, SO2.isa.
            d = size(tr);
            if ndims(tr) >= 2
                h =  all(d(1:2) == [4 4]);
                
                if h && nargin > 1 && ~isa(tr, 'sym')
                    h = SO3.isa( tr(1:3,1:3) );  % test the rotation submatrix
                    h = h && all(tr(4,:) == [0 0 0 1]);  % test the bottom row
                end
            else
                h = false;
            end
        end
        
        
    end
end
