clc;
clear;
close all;

% Create a Newton fractal for
% F1 =(x+y^2-4)
% F2 = 2x+y-2
%
% Subproblem solutions
% g = x+y^2 -4
% h = -2*y^2 + y +6

tol=1e-10;
npts=201; %npts=2001 gives good resolution
X=linspace(-2,2,npts);
Y=linspace(-2,2,npts);

[x,y]=meshgrid(X,X); % set up a grid of starting guesses for Newton's method

% define the three roots of f(z)
r1=[0  ;2];
r2 = [7/4; -1.5];
% r2=[3.2500000000000000000 , 2.7041634565979919698 ];
% r3=[3.2500000000000000000  -2.7041634565979919698  ];
% r4=[ 0.55570088181281412131  0.28038027781483410146  ];
% r5=[ 0.55570088181281412131 -0.28038027781483410146  ];
% r6=[ 3.2215084844679104104  2.1157099960468046642  ];
% r7=[ 3.2215084844679104104  -2.1157099960468046642  ];
% r8=[ 3.2215084844679104104 +2.1157099960468046642  ];
% r9=[ 3.2215084844679104104 -2.1157099960468046642  ];

%
% Newton with jacobian of non-preconditioned system
cp = zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k]=mynewt(x(i,j),y(i,j),tol);
        %         norm(res-r1)
        cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure(1);

hold on;
subplot(2,2,1)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
xlabel('x')
ylabel('y')
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
% colorbar
title('Original Jacobian','Fontsize', 14)

figure(2)

hold on;
subplot(2,2,1)
xt = get(gca, 'XTick');
xlabel('x')
ylabel('y')
set(gca, 'FontSize', 14)
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Original Jacobian','Fontsize', 14)
colorbar

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
        %         r1
        %         res
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure(1)

subplot(2,2,2)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
xlabel('x')
ylabel('y')
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('True MSPIN Jacobian','Fontsize', 14)
% colorbar
figure(2)
subplot(2,2,2)
xt = get(gca, 'XTick');
xlabel('x')
ylabel('y')
set(gca, 'FontSize', 14)
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('True MSPIN Jacobian','Fontsize', 14)
colorbar


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
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end

figure(1)
subplot(2,2,3)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
xlabel('x')
ylabel('y')
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Approx. MSPIN Jacobian M','Fontsize', 14)
% colorbar
figure(2)
subplot(2,2,3)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
xlabel('x')
ylabel('y')
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Approx. MSPIN Jacobian M','Fontsize', 14)
colorbar


%
% bad approximation to mspin jacobian


cp=zeros(size(x));
cit=zeros(size(x));
sing_mat = zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        
        [res, k,sing]=mynewt_mspinaprox2(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        sing_mat(i,j) = sing;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure(1)
subplot(2,2,4)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight

xlabel('x')
ylabel('y')
title('Approx. MSPIN Jacobian NM','Fontsize', 14)
hold off;
% colorbar
figure(2)
subplot(2,2,4)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Approx. MSPIN Jacobian NM','Fontsize',14)
colorbar
hold off;


[res_mspinex, k,path_mspinex]=mynewt_mspinexact_path(0,0,tol);
[res_newt, k,path_newt]=mynewt_path(0,0,tol);
[res_mspinap1, k,path_mspinap1]=mynewt_mspinaprox1_path(0,0,tol);
[res_mspinap2, k,path_mspinap2]=mynewt_mspinaprox2_path(0,0,tol);

for i = 1: size(path_mspinex,1)
    res_mspinex(i) = norm(f(path_mspinex(i,:)),2);
    res_mspinexpre(i) = norm(f_mspin(path_mspinex(i,:)),2);
end
for i = 1: size(path_newt,1)
    res_newt(i) = norm(f(path_newt(i,:)),2);
end
for i = 1: size(path_mspinap1,1)
    res_mspinap1(i) = norm(f(path_mspinap1(i,:)),2);
    res_mspinap1pre(i) = norm(f_mspin(path_mspinap1(i,:)),2);
end
for i = 1: size(path_mspinap2,1)
    res_mspinap2(i) = norm(f(path_mspinap2(i,:)),2);
    res_mspinap2pre(i) = norm(f_mspin(path_mspinap2(i,:)),2);
end

figure;
subplot(1,2,1)
hold on;
plot(1:length(res_mspinex),res_mspinex,'b:x')
plot(1:length(res_newt),res_newt,'r*')
plot(1:length(res_mspinap1),res_mspinap1,'k--o')
plot(1:length(res_mspinap2),res_mspinap2,'g--+')
xlabel('Iterate')
legend('Exact MSPIN','Newton','Approx. MSPIN M','Approx. MSPIN NM','Location','NorthEast')
ylabel('|| Residual||_{2}')
title('Residuals of the Original System')
hold off;
subplot(1,2,2)
hold on;
plot(1:length(res_mspinexpre),res_mspinexpre,'r:x')
plot(1:length(res_mspinap1pre),res_mspinap1pre,'k--o')
plot(1:length(res_mspinap2pre),res_mspinap2pre,'g--+')
xlabel('Iterate')
legend('Exact MSPIN','Approx. MSPIN M','Approx. MSPIN NM','Location','NorthEast')
ylabel('|| Residual||_{2}')
title('MSPIN Residuals')
hold off
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


% function for good mspin approx
function [xnew,k,path]=mynewt_mspinaprox1_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
path = [];
while (norm(fxold,Inf)>tol)
    
    path(k+1,1:2)= xold';
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


% function for good mspin approx
function [xnew,k,path]=mynewt_mspinaprox2_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
path = [];
while (norm(fxold,Inf)>tol)
    
    path(k+1,1:2)= xold';
    xnew=xold-(Japrox2(xold))\fxold;
    xold=xnew;
    fxold=f_mspin(xold);
    k=k+1;
    if k > 50
        break;
    end
end
xnew=xold;
end

% Path function for original newton
function [xnew,k,path]=mynewt_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f(xold);
k=0;
path = [];
while (norm(fxold,Inf)>tol)
    
    path(k+1,1:2)= xold';
    xnew=xold-(Jex(xold))\fxold;
    xold=xnew;
    fxold=f(xold);
    k=k+1;
    if k > 50
        break;
    end
end
xnew=xold;
end

% function for exact mspin jacobian
function [xnew,k,path]=mynewt_mspinexact_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
path = [];
while (norm(fxold,Inf)>tol)
    
    path(k+1,1:2)= xold';
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
function fval=f(b)
x=b(1); y = b(2);
% evaluate function
fval=[(x+y^2 -4); (2*x + y - 2)];
end


% function for goold approximate mspin jacobian
function [xnew,k,sing]=mynewt_mspinaprox2(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
sing = 0;
while (norm(fxold,Inf)>tol)
    Jap = Japrox2(xold);
    if det(Jap)==0
        sing = 1;
    end
    xnew=xold-(Jap)\fxold;
    xold=xnew;
    fxold=f_mspin(xold);
    k=k+1;
    if k > 50
        break;
    end
end
xnew=xold;
end

function fval=f_mspin(x)
% evaluate function
fval=[(x(1)+x(2)^2 -4); (-2*x(2)^2+6+x(2))];
end


function J_exact = Jex(x)

J_exact= [1, 2*x(2); 2 1];
end

function J_mspin = Jmspin(x)
J_mspin= [1 2*x(2);  0 (-4*x(2) + 1)] ;
end

function J_mspin_aprox = Japrox(x)
gh=f_mspin(x);
% display('gh is\n')
% gh
% caculate [p;v]=[u-g;v] (x_soln=[u,v])
pv=x-[gh(1);0];
J = Jex(pv); %Jex(x)
L = tril(J);
J_mspin_aprox = L\J;
end

function J_mspin_aprox = Japrox2(x)

J = Jex(x);
L = tril(J);
J_mspin_aprox = L\J;
end
