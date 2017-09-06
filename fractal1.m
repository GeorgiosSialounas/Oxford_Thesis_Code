clc;
clear;
close all;

% Create a Newton fractal for
% F1 =(x-y^3+1)^3-y^3
% F2 = 2x+3y-5
%
% Subproblem solutions
% g = x-y^3+1 -y
% h = 2/3*y^3+5/3*y-7/3

tol=1e-10;
npts=201; %npts=2001 gives good resolution
X=linspace(-1,3,npts);
Y=linspace(-1,3,npts);

[x,y]=meshgrid(X,X); % set up a grid of starting guesses for Newton's method

% define the three roots of f(z)
r1=[1  ;1];
% r2=[3.2500000000000000000 , 2.7041634565979919698 ];
% r3=[3.2500000000000000000  -2.7041634565979919698  ];
% r4=[ 0.55570088181281412131  0.28038027781483410146  ];
% r5=[ 0.55570088181281412131 -0.28038027781483410146  ];
% r6=[ 3.2215084844679104104  2.1157099960468046642  ];
% r7=[ 3.2215084844679104104  -2.1157099960468046642  ];
% r8=[ 3.2215084844679104104 +2.1157099960468046642  ];
% r9=[ 3.2215084844679104104 -2.1157099960468046642  ];

% 
% % Newton with jacobian of non-preconditioned system
% cp=zeros(size(x));
% cit=zeros(size(x));
% % run Newton's method for each starting guess
% for i=1:npts
%     %i
%     for j=1:npts
%         [res, k]=mynewt(x(i,j),y(i,j),tol);
%         norm(res-r1)
%         cit(i,j) = k;
%         if norm(res-r1)<tol       % if there is convergence to root 1
%             cp(i,j)=1;            % then record this
%             cmm(i,j,:)=[1,0,0];   % record colour to plot as red
%         
%         else                      % if there is no convergence
%             cmm(i,j,:)=[1,1,1];   % record colour to plot as white
%         end
%     end
% end
% figure;
% surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
% title('Exact: root')
% 
% figure;
% surf(x,y,cit), view(2), shading interp, axis equal tight
% title('Exact: iteration')



% Newton with exact MSPIN jacobian
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k]=mynewt_mspinexact(x(i,j),y(i,j),tol);
        cit(i,j) = k;
%         norm(res-r1)
r1
res
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
        
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure;
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
title('MSPIN Exact: root')

figure;
surf(x,y,cit), view(2), shading interp, axis equal tight
title('MSPIN Exact: iteration')

% Newton with approximate MSPIN jacobian
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k]=mynewt_mspinaprox(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
        
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure;
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
title('MSPIN approx: root')

figure;
surf(x,y,cit), view(2), shading interp, axis equal tight
title('MSPIN approx: iteration')


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
fval=[((x(1)-3*x(2)^3 +1)^3 -x(2)^3); (2*x(1)+3*x(2)-5)];
end

function fval=f_mspin(x)
% evaluate function
fval=[((x(1)-3*x(2)^3 +1) -x(2));( (2/3)*x(2)^3 + (5/3)*x(2)-7/3)];
end


function J_exact = Jex(x)
J_exact= [3*(x(1)-x(2)^3+1)^2 , ((3*(x(1)-x(2)^3+1)^2)*(-3*x(2)^2)-3*x(2)^2); 2 3];
end

function J_mspin = Jmspin(x)
J_mspin= [1 (-3*x(2)^2 -1);  0 (2*x(2)^2+5/3)] ;
end

function J_mspin_aprox = Japrox(x)
J = Jex(x);
L = tril(J);
J_mspin_aprox = L\J;
end

