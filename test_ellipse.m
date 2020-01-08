mu = [0;0]
Sigma = [2 3;
         3 10];
     
level = 1;
npoints = 100;


[ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )



[V,D] = eig(Sigma);
ang = linspace(0,2*pi,npoints);
circle = [cos(ang); sin(ang)];
ell = V*level*sqrt(D)*circle + mu;



figure;
hold on; grid on; axis equal
plot(xy(1,:),xy(2,:))
plot(ell(1,:),ell(2,:),'--')

quiver(mu(1),mu(2),sqrt(D(1,1))*V(1,1),sqrt(D(1,1))*V(2,1))
quiver(mu(1),mu(2),sqrt(D(2,2))*V(1,2),sqrt(D(2,2))*V(2,2))