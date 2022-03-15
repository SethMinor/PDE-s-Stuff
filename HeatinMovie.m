clear, clc, clf;

% rectangular spatial domain parameters
L = 10; % length
H = 5; % height

% temporal parameters 
k = 0.02; % diffusivity
tf = 100; % final time
num_frames = 150; % number of movie frames
time = linspace(20,tf,num_frames); % time span

% other parameters
p = 5.1;
q = p;

% spatial eigenfunctions 
phi =@(n,m,x,y) sin(((n*pi)/L).*x).*sin(((m*pi)/H).*y); % \phi_{n,m}

% temporal part of solution 
h =@(n,m,t) exp(-k*(((n*pi)/L)^2+((m*pi)/H)^2)*t);

% creating the spatial mesh 
nx = 60; % resolution of the mesh
ny = 60;

rows = 5; % n's
columns = 5; % m's

x = linspace(0,L,nx);
y = linspace(0,H,nx);
[X,Y] = meshgrid(x,y);

% defining the coefficients (given L,H,p,q)
cn =@(n) (4*L*p^2*sin(n*pi/2)*sin((n*pi)/(2*p)))/(n*pi*(4*p^2-n^2)); 
cm =@(m) (4*H*q^2*sin(m*pi/2)*sin((m*pi)/(2*q)))/(m*pi*(4*q^2-m^2));
Anm =@(n,m) ((pi^2)/(L*H))*cn(n)*cm(m);

%% create the movie by plotting the sum of the spatial modes
figgy = figure;

axis = gca; 
colormap(axis,hot(30))
axis.NextPlot = 'replaceChildren';
grid on
bagOframes(num_frames) = struct('cdata',[],'colormap',[]);

%V = VideoWriter('Diffusing_Heat.avi');
%open(V);

% increment over time
j = 1;
for t=time

    U = 0; % initialize the solution

    % increment over space
    for n = 1:rows % n's
        for m = 1:columns % m's
            A_nm = Anm(n,m); % compute coefficients
            
            phi_nm = phi(n,m,X,Y); % spatial modes
            h_nm = h(n,m,t); % temporal solution

            U = U + phi_nm*h_nm; % add another (n,m) mode
        end
    end

    % plot solution at time t
    surf(X,Y,U)
    zlim([-1 3])
    subtitle("Time: " + t,'interpreter','latex')
    xlabel('$x$ (for $L=10$)','interpreter','latex')
    ylabel('$y$ (for $H=5$)','interpreter','latex')
    zlabel('$u = u(x,y,t)$','interpreter','latex')
    drawnow
    pause(0.01)
    
    frame = getframe(axis);
    bagOframes(j) = frame;
    
    %writeVideo(V,frame);

    j = j + 1;

end

%close(V);

% play the movie
movie(bagOframes,1,30);


