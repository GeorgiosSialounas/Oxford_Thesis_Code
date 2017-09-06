clc;
clear;
close all;

% G=F1(x1,x2)=(x1-x2^3+1)^3-x2^3
% H=F2(x1,x2)=x1+2x2-3
% The exact solution is x_exact=[1,1]';

% g=x1-x2^3-x2+1
% h=0.5*x2^3+1.5*x2-2

tol=1e-3;
npts=201; %npts=2001 gives good resolution
X=linspace(-2,4,npts);
[x,y]=meshgrid(X,X); % set up a grid of starting guesses for Newton's method

% define the two roots of f(z)
r1=[1 1];

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
%         elseif norm(res-r2)<tol   % if there is convergence to root 2
%             cp(i,j)=2;            % then record this
%             cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
% %             disp('root 2')
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
%         elseif norm(res-r2)<tol   % if there is convergence to root 2
%             cp(i,j)=2;            % then record this
%             cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
% %             disp('root 2')
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
%         norm(res-r1)
%         disp('here')
%         i
%         j
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
%             disp ('root 1')
%         elseif norm(res-r2)<tol   % if there is convergence to root 2
%             cp(i,j)=2;            % then record this
%             cmm(i,j,:)=[0,1,0];   % record colour to plot as green            
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
    if k > 50
        break;
    end
end
xnew=xold;
end

function fval=f(x)
% evaluate function

fval=[((x(1)-x(2)^3+1)^3-x(2)^3); (x(1)+2*x(2)-3)];
end

function fval=f_mspin(x)
% evaluate function
fval=[(x(1)-x(2)-x(2)^3+1); (0.5*x(2)^3+1.5*x(2)-2)];
end



function J_exact = Jex(x)
J_exact= [3*(x(1)-x(2)^3+1)^2, -9*x(2)^2*(x(1)-x(2)^3+1)^2-3*x(2)^2;1,2];
end

function J_mspin = Jmspin(b)
x = b(1); y = b(2);
J_mspin= [1 (-3*y^2-1); 0 (1.5*y^2+1.5)];
end

function J_mspin_aprox = Japrox(x)
% J_exact= Jex(x);
% L = tril(J_exact);
% J_mspin_aprox = L\J_exact;
gh=f_mspin(x);
     % display('gh is\n')
     % gh
    % caculate [p;v]=[u-g;v] (x_soln=[u,v])
pv=x-[gh(1);0];
J = Jex(pv); %Jex(x)
L = tril(J);
J_mspin_aprox = L\J;
end

