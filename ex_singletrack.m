d = @(z)deal(0,0);
sigmaw = 0;

x0 = [-2.5;0;3;0;pi/2;0;0;0;0;0];
u0 = [ 0 1 0 0 0]';
obj = MotionModelGP_SingleTrack(d,sigmaw)

obj.f(x0,u0)