%TWIST SE(3) Twist class
%
% A Twist class holds the parameters of a twist, a representation of a
% rigid body displacement in SE(3).
%
% Methods::
%  S             twist vector (1x6)
%  se            twist as (augmented) skew-symmetric matrix (4x4)
%  T             convert to homogeneous transformation (4x4)
%  R             convert rotational part to matrix (23x3)
%  exp           synonym for T
%  ad            logarithm of adjoint
%  pitch         pitch of the screw
%  pole          a point on the line of the screw
%  prod          product of a vector of Twists
%  theta         rotation about the screw
%  line          Plucker line object representing line of the screw
%  display       print the Twist parameters in human readable form
%  char          convert to string
%
% Conversion methods::
%  SE            convert to SE3 object
%  double        convert to real vector
%
% Overloaded operators::
%  *             compose two Twists
%  *             multiply Twist by a scalar
%
% Properties (read only)::
%  v             moment part of twist (3x1)
%  w             direction part of twist (3x1)
%
% References::
% - "Mechanics, planning and control"
%   Park & Lynch, Cambridge, 2016.
%
% See also trexp, trexp2, trlog.

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

classdef Twist
    properties (SetAccess = protected)
        v  %axis direction (column vector)
        w  %moment (column vector)
    end
    
    methods
        function tw = Twist(T, varargin)
            %Twist.Twist Create Twist object
            %
            % TW = Twist(T) is a Twist object representing the SE(2) or SE(3)
            % homogeneous transformation matrix T (3x3 or 4x4).
            %
            % TW = Twist(V) is a twist object where the vector is specified directly.
            %            %
            % TW = Twist('R', A, Q) is a Twist object representing rotation about the
            % axis of direction A (3x1) and passing through the point Q (3x1).
            %
            % TW = Twist('R', A, Q, P) as above but with a pitch of P (distance/angle).
            %
            % TW = Twist('T', A) is a Twist object representing translation in the
            % direction of A (3x1).
            %
            % Notes::
            %  The argument 'P' for prismatic is synonymous with 'T'.
            
            if ischar(T) || isstring(T)
                % 'P', dir
                % 'R', dir, point 3D
                % 'R', point   2D
                switch upper(char(T))
                    case 'R'
                        
                        dir = varargin{1};
                        if length(dir) < 3
                            error('SMTB:Twist:badarg', 'For 2d case can only specify position');
                        end
                        point = varargin{2};
                        
                        w = unit(dir(:))';
                        
                        v = -cross(w, point(:)');
                        if nargin >= 4
                            pitch = varargin{3};
                            v = v + pitch * w;
                        end
                        
                        
                    case {'P', 'T'}
                        dir = varargin{1};
                        
                        w = [0 0 0];
                        v = unit(dir(:)');
                end
                
                if ~isa(v, 'sym')
                    v(abs(v)<eps) = 0;
                end
                if ~isa(w, 'sym')
                    w(abs(w)<eps) = 0;
                end
                tw.v = v;
                tw.w = w;
            elseif size(T,1) == size(T,2)
                % it's a square matrix
                if T(end,end) == 1
                    % its a homogeneous matrix, take the logarithm
                    S = logm(T);  % use closed form for SE(3)

                    [skw,v] = tform2rotm(S);
                    tw.v = v';
                    tw.w = skew2vec(skw);
                else
                    % it's an augmented skew matrix, unpack it
                    [skw,v] = tform2rotm(T);
                    tw.v = v';
                    tw.w = skew2vec(skw);
                end
            elseif isvector(T)
                % its a row vector form of twist, unpack it
                if length(T) == 6
                    tw.v = T(1:3)'; tw.w = T(4:6)';
                else
                    error('SMTB:Twist:badarg', '3 or 6 element vector expected');
                end
            end
        end

        function printline(obj, varargin)
            trprint(obj.exp(1), varargin{:});
        end
        
        function Su = unit(S)
            %Twist.unit Return a unit twist
            %
            % TW.unit() is a Twist object representing a unit aligned with the Twist
            % TW.
            if abs(S.w) > 10*eps
                % rotational twist
                Su = Twist( double(S) / norm(S.w) );
            else
                % prismatic twist
                Su = Twist( [unit(S.v); 0; 0; 0] );
            end
        end
        
        function x = S(tw)
            %Twist.S Return the twist vector
            %
            % TW.S is the twist vector in se(3) as a vector (6x1).
            %
            % Notes::
            % - Sometimes referred to as the twist coordinate vector.
            x = [tw.v tw.w];
        end
        
        function x = double(tw)
            %Twist.double Return the twist vector
            %
            % double(TW) is the twist vector in se(3) as a vector (6x1). 
            % If TW is a vector (1xN) of Twists the result is a matrix (6xN) with
            % one column per twist.
            %
            % Notes::
            % - Sometimes referred to as the twist coordinate vector.
            x = [tw.v tw.w];
        end
        
        function x = se3(tw)
            %Twist.se Return the twist matrix
            %
            % TW.se is the twist matrix in sse(3) which is an augmented
            % skew-symmetric matrix (4x4).
            %
            x = vec2skewa(tw.S);
        end
        
        
        function c = mtimes(a, b)
            %Twist.mtimes Multiply twist by twist or scalar
            %
            % TW1 * TW2 is a new Twist representing the composition of twists TW1 and
            % TW2.
            %
            % TW * T is an SE3 that is the composition of the twist TW and the
            % homogeneous transformation object T.
            %
            % TW * S with its twist coordinates scaled by scalar S.
            %
            % TW * T compounds a twist with an SE3 transformation
            %
            
            if isa(a, 'Twist')
                if isa(b, 'Twist')
                    % twist composition
                    c = Twist( a.exp * b.exp);
                elseif isreal(b)
                    c = Twist(a.S * b);
                elseif length(a.v) == 3 && ishomog(b)
                    % compose a twist with SE3, result is an SE3
                    c = SE3(a.T * double(b));
                elseif isa(b, 'SpatialVelocity')
                    c = SpatialVelocity(a.Ad * b.vw);
                elseif isa(b, 'SpatialAcceleration')
                    c = SpatialAcceleration(a.Ad * b.vw);
                elseif isa(b, 'SpatialForce')
                    c = SpatialForce(a.Ad' * b.vw);
                else
                    error('SMTB:Twist', 'twist * SEn, operands don''t conform');
                end
            elseif isreal(a) && isa(b, 'Twist')
                c = Twist(a * b.S);
            elseif isa(a, 'Twist') && isreal(b)
                c = Twist(a.S * b);
            else
                error('SMTB:Twist: incorrect operand types for * operator')
            end
        end
        
        function x = mrdivide(a, b)
            x = Twist(a.S / b);
        end
            
        function x = SE3(tw, varargin)
            %Twist.exp Convert twist to homogeneous transformation
            %
            % TW.SE3 is the homogeneous transformation equivalent to the twist (SE2 or SE3).
            %
            % TW.SE3(THETA) as above but with a rotation of THETA about the twist.
            %
            % Notes::
            % - For the second form the twist must, if rotational, have a unit rotational component.
            %
            % See also Twist.T, trexp, trexp2.
            opt.deg = false;
            [opt,args] = tb_optparse(opt, varargin);
            
            if opt.deg && all(tw.w == 0)
                warning('Twist: using degree mode for a prismatic twist');
            end
            
            if ~isempty(args)
                theta = args{1};
                
                if opt.deg
                    theta = theta * pi/180;
                end
            else
                theta = 1;
            end
            
            ntheta = length(theta);
            assert(length(tw) == ntheta || length(tw) == 1, 'Twist:exp:badarg', 'length of twist vector must be 1 or length of theta vector')
            
            x(ntheta) = SE3;
            if length(tw) == ntheta
                for i=1:ntheta
                    x(i) = expm( vec2skewa(tw(i).S * theta(i)) );
                end
            else
                for i=1:ntheta
                    x(i) = expm( vec2skewa(tw.S * theta(i)) );
                end
            end
        end
        
        function x = ad(tw)
            %Twist.ad Logarithm of adjoint
            %
            % TW.ad is the logarithm of the adjoint matrix of the corresponding
            % homogeneous transformation.
            %
            % See also SE3.Ad.
            x = [ vec2skew(tw.w) vec2skew(tw.v); zeros(3,3) vec2skew(tw.w) ];
        end
        
        function x = Ad(tw)
            %Twist.Ad Adjoint
            %
            % TW.Ad is the adjoint matrix of the corresponding
            % homogeneous transformation.
            %
            % See also SE3.Ad.
            x = tw.SE.Ad;
        end
        
        
        function out = SE(tw)
            %Twist.SE Convert twist to SE3 object
            %
            % TW.SE is an SE3 object representing the homogeneous transformation equivalent to the twist.
            %
            % See also Twist.T, SE3.
            if length(tw.v) == 2
                out = SE2( tw.T );
            else
                out = SE3( tw.T );
            end
        end
        
        function x = exp(tw, varargin)
            %Twist.T Convert twist to homogeneous transformation
            %
            % TW.T is the homogeneous transformation equivalent to the twist (4x4).
            %
            % TW.T(THETA) as above but with a rotation of THETA about the twist.
            %
            % Notes::
            % - For the second form the twist must, if rotational, have a unit rotational component.
            %
            % See also Twist.exp, trexp, trexp2, trinterp, trinterp2.
%             x = double( tw.exp(varargin{:}) );
            x = expm(tw.se3);
        end
        
        
        function p = pitch(tw)
            %Twist.pitch Pitch of the twist
            %
            % TW.pitch is the pitch of the Twist as a scalar in units of distance per radian.

            p = tw.w * tw.v';
        end
        
        function L = line(tw)
            %Twist.line Line of twist axis in Plucker form
            %
            % TW.line is a Plucker object representing the line of the twist axis.
            %
            % Notes::
            % - For 3D case only.
            %
            % See also Plucker.
            
            % V = -tw.v - tw.pitch * tw.w;
            for i=1:length(tw)
                L(i) = Plucker([ -tw(i).v - tw(i).pitch * tw(i).w tw(i).w] ); %#ok<AGROW>
            end
        end
        
        function out = prod(obj)
            %Twist.prod Compound array of twists
            %
            % TW.prod is a twist representing the product (composition) of the
            % successive elements of TW (1xN), an array of Twists.
            %
            %
            % See also RTBPose.prod, Twist.mtimes.
            out = obj(1);
            
            for i=2:length(obj)
                out = out * obj(i);
            end
        end
        
        function p = pole(tw)
            %Twist.pole Point on the twist axis
            %
            % TW.pole is a point on the twist axis (3x1).
            %
            % Notes::
            % - For pure translation this point is at infinity.

            p = cross(tw.w, tw.v) / tw.theta();
        end
        
        function th = theta(tw)
            %Twist.theta Twist rotation
            %
            % TW.theta is the rotation (1x1) about the twist axis in radians.
            %
            
            th = norm(tw.w);
        end
        
        
        function s = char(tw)
            %Twist.char Convert to string
            %
            % s = TW.char() is a string showing Twist parameters in a compact single line format.
            % If TW is a vector of Twist objects return a string with one line per Twist.
            %
            % See also Twist.display.
            s = '';
            for i=1:length(tw)
                
                ps = '( ';
                ps = [ ps, sprintf('%0.5g  ', tw(i).v) ]; %#ok<AGROW>
                ps = [ ps(1:end-2), '; '];
                ps = [ ps, sprintf('%0.5g  ', tw(i).w) ]; %#ok<AGROW>
                ps = [ ps(1:end-2), ' )'];
                if isempty(s)
                    s = ps;
                else
                    s = char(s, ps);
                end
            end
            
            
        end
        
        function display(tw) %#ok<DISPLAY>
            %Twist.display Display parameters
            %
            % L.display() displays the twist parameters in compact single line format.  If L is a
            % vector of Twist objects displays one line per element.
            %
            % Notes::
            % - This method is invoked implicitly at the command line when the result
            %   of an expression is a Twist object and the command has no trailing
            %   semicolon.
            %
            % See also Twist.char.
            loose = strcmp( get(0, 'FormatSpacing'), 'loose'); %#ok<GETFSP>
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(tw) );
        end % display()
        
    end

    methods(Static)
        function tw = UnitRevolute(varargin)
            tw = Twist('R', varargin{:});
        end

        function tw = UnitPrismatic(varargin)
            tw = Twist('P', varargin{:});
        end
    end
end
