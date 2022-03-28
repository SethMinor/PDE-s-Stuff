% Implementing the Crank-Nicolson FDS as applied
% to the heat equation, using LU decomp
clear, clc, clf;

% time and space discretization
h = 1/40;
%lambda = 1;
%k = lambda*h;
mu = 10;
k = mu*h^2;
lambda = k/h;


% initial and final times
t0 = 0;
tf = 1/2;

% time step stuff
time = t0:k:tf;
N = length(time);

% space step stuff
L = -1;
R = 1;
x = L:h:R;
M = length(x);

% initialize numerical solution vectors
V = zeros(M,1);

% inish condishes at time t=t0 (or n=0)
for j = 1:M
    if (abs(x(j)) < 1/2)
        V(j) = 1;
    elseif (abs(x(j)) == 1/2)
        V(j) = 1/2;
    else 
        V(j) = 0;
    end
end

% more constants 
B = 1; % diffusion coeff
alpha = B/(2*h^2);
beta = 1/k + B/h^2; 

% analytical solution
cap = 1000; % truncation of infinite sum
Uexact =@(X,T) exact(X,T,cap);



% create L and U matrices
A = zeros(M);
for i=1:M
    for j=1:M
        if (((i+1) == j) && (i ~= 1))
            A(i,j) = -alpha;
        elseif ((i == (j+1)) && (i ~= M))
            A(i,j) = -alpha;
        elseif (i == j)
            A(i,j) = beta;
        else
            A(i,j) = 0;
        end
    end
end
A(1,1) = 1;
A(M,M) = 1;

% let the proprietary software do its thing
[L,U] = lu(A);
Linv = inv(L);
Uinv = inv(U);
updateMatrix = Uinv*Linv;

% initialize the time counter
t = t0;

figure (1)

plot(x,V,'-x')
xlabel('$x$','interpreter','latex')
ylabel('$u = u(x,t)$','interpreter','latex')
title("$u = u(x,t)$ at $t=$ "+t+...
    ", for $h=$ "+h+" and $\lambda=$"+lambda+...
    ", $\mu=$ "+mu,...
    'interpreter','latex')
grid on
hold on
plot(x,Uexact(x,t))
legend('Numerical solution','Exact solution',...
    'interpreter','latex')
ylim([0 1.2])
pause(0.5)
hold off

% main time loop
for n=1:N-1
    % update time counter 
    t = t + k;

    % make b vector
    b = makeb(V,M,alpha,beta,(n-1),k,cap);
    
    % find V_new using L,U and b
    V = updateMatrix*b;

    % plot stuff
    plot(x,V,'-x')
    xlabel('$x$','interpreter','latex')
    ylabel('$u = u(x,t)$','interpreter','latex')
    title("$u = u(x,t)$ at $t=$ "+t+...
    ", for $h=$ "+h+" and $\lambda=$"+lambda+...
    ", $\mu=$ "+mu,...
    'interpreter','latex')
    ylim([0 1.2])
    grid on
    hold on
    %subplot(1,2,2)
    plot(x,Uexact(x,t))
    legend('Numerical solution','Exact solution',...
    'interpreter','latex')
    pause(0.01)
    hold off

end
%%
function u = exact(x,t,cap)
% This function evaluates the exact solution at (x,t) by approximating
% it using a partial sum running from l=0 to l=cap.
    u = 0;
        for l=0:cap
            u = u +...
                ((-1)^l*(cos(pi*(2*l+1).*x))./(pi*(2*l+1)))...
                *exp(-pi^2*(2*l+1).^2.*t);
        end
    u = 1/2 + 2*u;

end

%%
function bn = makeb(Vn,M,alpha,beta,n,k,cap)
% This returns the b^n vector that is a function of a corresponding
% V^n vector
    bn = zeros(M,1);
    % boundary conditions
    bn(1) = exact(-1,(n+1)*k,cap);
    bn(M) = exact(1,(n+1)*k,cap);
    for m=2:M-1
        bn(m) = alpha*Vn(m+1) + (beta - 4*alpha)*Vn(m) + alpha*Vn(m-1);
    end
end