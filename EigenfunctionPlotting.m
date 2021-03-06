clear, clc, clf;

% rectangular spatial domain parameters
L = 10; % length
H = 5; % height

% eigenvalues
mu =@(n) ((n*pi)/L).^2; % \mu_{n}
lambda =@(n,m) ((n*pi)/L).^2 + ((m*pi)/H).^2; % \lambda_{n,m}

% spatial eigenfunctions 
phi =@(n,m,x,y) sin(((n*pi)/L).*x).*sin(((m*pi)/H).*y); % \phi_{n,m}

% creating the spatial mesh 
nx = 30; % resolution of the mesh
ny = 30;

rows = 5; % subplot stuff 
columns = 5; 

x = linspace(0,L,nx);
y = linspace(0,H,nx);
[X,Y] = meshgrid(x,y);

% plot the spatial modes 
figure (1)

for i = 1:rows
    for j = 1:columns
        subplot(rows, columns, j+(i-1)*columns)
        Z = phi(i,j,X,Y);
        surfc(X,Y,Z)
        title("$(n,m)= ($" + i + "," + j + ")",'interpreter','latex')
    end
end

sgtitle('Spatial eigenfunctions $\phi_{nm}$ for $n,m  = 1, \dots, 5$',...
    'interpreter','latex')