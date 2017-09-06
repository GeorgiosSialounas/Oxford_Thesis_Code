clc;
clear;
close all;

% Create a Newton fractal for
% F1 = (x+y^2)-4
% F2 = 2x+y-2
%
% Subproblem solutions
% g = x+y^2-4
% h = -2y^2+y+6

tol=1e-3;
npts=501; %npts=2001 gives good resolution
X=linspace(-4,4,npts);
[x,y]=meshgrid(X,X); % set up a grid of starting guesses for Newton's method

% define the two roots of f(z)
r1=[0;  2];
r2=[3; -1 ];

cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
%             disp ('root 1')
        elseif norm(res-r2)<tol   % if there is convergence to root 2
            cp(i,j)=2;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
%             disp('root 2')
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end

figure(1)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
title('Root: Exact')

figure(2)
surf(x,y,cit), view(2), shading interp, axis equal tight
colorbar
title('Iteration exact')


% Exact MSPIN
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt_mspinexact(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
%             disp ('root 1')
        elseif norm(res-r2)<tol   % if there is convergence to root 2
            cp(i,j)=2;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
%             disp('root 2')
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end

figure(3)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
title('Root: MSPIN Exact')

figure(4)
surf(x,y,cit), view(2), shading interp, axis equal tight
colorbar
title('Iteration: MSPIN exact')

% Aprox MSPIN
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt_mspinaprox(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
%             disp ('root 1')
        elseif norm(res-r2)<tol   % if there is convergence to root 2
            cp(i,j)=2;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
%             disp('root 2')
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end

figure(5)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
title('Root: MSPIN aprox')

figure(6)
surf(x,y,cit), view(2), shading interp, axis equal tight
colorbar
title('Iteration: MSPIN aprox')



function [xnew,k]=mynewt(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f(xold);
k=0;
while (norm(fxold,Inf)>tol)
    xnew=xold-(Jex(xold))\fxold;
    xold=xnew;
    fxold=f(xold);
    k=k+1;
    if k >50
        break;
    end
end
xnew=xold;
end

% function for exact mspin jacobian
function [xnew,k]=mynewt_mspinexact(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
while (norm(fxold,Inf)>tol)
    xnew=xold-(Jmspin(xold))\fxold;
    xold=xnew;
    fxold=f_mspin(xold);
    k=k+1;
    if k > 50
        break;
    end
end
xnew=xold;
end

% function for approximate mspin jacobian
function [xnew,k]=mynewt_mspinaprox(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
while (norm(fxold,Inf)>tol)
    xnew=xold-(Japrox(xold))\fxold;
    xold=xnew;
    fxold=f_mspin(xold);
    k=k+1;
end
xnew=xold;
end

function fval=f(b)
x = b(1); y = b(2);
% evaluate function
fval=[(x + y^2 -4); (2*x+ y -2)];
end

function fval=f_mspin(b)
x = b(1); y = b(2);
% evaluate function
fval=[(x + y^2 -4); (-y^2 + y -2)];
end

function Jval=J(x)
% evaluate Jacobian
t1=3*(x(1)^2-x(2)^2);
t2=6*x(1)*x(2);
Jval=[t1 -t2; t2 t1];
end

function J_exact = Jex(b)
x = b(1); y = b(2);
J_exact= [1 2*y; 1 1];
end

function J_mspin = Jmspin(b)
x = b(1); y = b(2);
J_mspin= [1 2*y; 0 -2*y+1];
end

function J_mspin_aprox = Japrox(b)
% x = b(1); y = b(2);
% L = [1 0; 1 1];
% Jm = [1 2*y; 1 1];
% J_mspin_aprox = L\Jm;
gh=f_mspin(x);
     % display('gh is\n')
     % gh
    % caculate [p;v]=[u-g;v] (x_soln=[u,v])
pv=x-[gh(1);0];
J = Jex(pv); %Jex(x)
L = tril(J);
J_mspin_aprox = L\J;s

end

