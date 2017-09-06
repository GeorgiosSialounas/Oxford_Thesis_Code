clc;
clear;
close all;

% Create a Newton fractal for
% F1 = ((x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -9*y+32)^3- y^3); 
% F2 = (x+y + 2)
%
% Subproblem solutions
% g = (x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -10*y+32); 
% h = (y^5 - 4*y^4 - 12*y^3 +34*y^2 +11*y -30)

tol=1e-10;
npts=201; %npts=2001 gives good resolution
X=linspace(-9,3,npts);
Y=linspace(-5,7,npts);

[x,y]=meshgrid(X,Y); % set up a grid of starting guesses for Newton's method

% define the three roots of f(z)
r1=[-7  ; 5];
r2 =[-4; 2];
r3 = [-3; 1];
r4 = [-1;-1];
r5 = [1;-3];
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
cit= 50*zeros(size(x)); % 
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k]=mynewt(x(i,j),y(i,j),tol);
        %         norm(res-r1)
%         cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
          elseif norm(res-r3)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0,1];   % record colour to plot as red
        elseif norm(res-r4)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0.5,0.5,0];   % record colour to plot as red
          elseif norm(res-r5)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0.5,0.5];   % record colour to plot as red    
            
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
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
% colorbar
title('Original Jacobian','Fontsize', 14)

figure(2)

hold on;
subplot(2,2,1)
xt = get(gca, 'XTick');
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
 %         cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
          elseif norm(res-r3)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0,1];   % record colour to plot as red
        elseif norm(res-r4)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0.5,0.5,0];   % record colour to plot as red
          elseif norm(res-r5)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0.5,0.5];   % record colour to plot as red    
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end
figure(1)

subplot(2,2,2)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('True MSPIN Jacobian','Fontsize', 14)
% colorbar
figure(2)
subplot(2,2,2)
xt = get(gca, 'XTick');
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
%         cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
          elseif norm(res-r3)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0,1];   % record colour to plot as red
        elseif norm(res-r4)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0.5,0.5,0];   % record colour to plot as red
          elseif norm(res-r5)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0.5,0.5];   % record colour to plot as red    
            
        else                      % if there is no convergence
            cmm(i,j,:)=[1,1,1];   % record colour to plot as white
        end
    end
end

figure(1)
subplot(2,2,3)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Approx. MSPIN Jacobian M','Fontsize', 14)
% colorbar
figure(2)
subplot(2,2,3)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
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
 %         cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
        elseif norm(res-r2)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,1,0];   % record colour to plot as red
          elseif norm(res-r3)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0,1];   % record colour to plot as red
        elseif norm(res-r4)<tol       % if there is convergence to root 1
            cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0.5,0.5,0];   % record colour to plot as red
          elseif norm(res-r5)<tol       % if there is convergence to root 1
              cit(i,j) = k;
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[0,0.5,0.5];   % record colour to plot as red    
            
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



[res_mspinex, k,path_mspinex]=mynewt_mspinexact_path(-0.5,-0.5,tol);
[res_newt, k,path_newt]=mynewt_path(-0.5,-0.5,tol);
[res_mspinap1, k,path_mspinap1]=mynewt_mspinaprox1_path(-0.5,-0.5,tol);
[res_mspinap2, k,path_mspinap2]=mynewt_mspinaprox2_path(-0.5,-0.5,tol);

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
plot(1:length(res_mspinex),res_mspinex,'b--*')
plot(1:length(res_newt),res_newt,'r:*')
plot(1:length(res_mspinap1),res_mspinap1,'k--o')
plot(1:length(res_mspinap2),res_mspinap2,'m--+')
xlabel('Iterate')
legend('Exact MSPIN','Newton','Approx. MSPIN M','Approx. MSPIN NM','Location','NorthEast')
ylabel('|| Residual||_{2}')
title('Original System Residuals','Fontsize',16)
hold off;
subplot(1,2,2)
hold on;
plot(1:length(res_mspinexpre),res_mspinexpre,'b--*')
plot(1:length(res_mspinap1pre),res_mspinap1pre,'k--o')
plot(1:length(res_mspinap2pre),res_mspinap2pre,'m--+')
xlabel('Iterate')
legend('Exact MSPIN','Approx. MSPIN M','Approx. MSPIN NM','Location','NorthEast')
ylabel('|| Residual||_{2}')
title('MSPIN Residuals','Fontsize',16)
hold off

F = zeros(size(x));
FMSPIN = zeros(size(x));
for  i = 1 : npts
    for j = 1 : npts
        F(i,j) = log(norm(f([x(i,j),y(i,j)]),2)+1);
        FMSPIN(i,j) = log(norm(f_mspin([x(i,j),y(i,j)]),2)+1);
   
    end
end
[~,~,pathap1]=mynewt_mspinaprox1_path(0,0, tol);
[~,~,pathap2]=mynewt_mspinaprox2_path(0,0, tol);

figure;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
subplot(1,2,1)
contour(x,y,F)
xlabel('x')
ylabel('y')
colorbar
title('log(||F(x) +1 ||)','Fontsize',16)
subplot(1,2,2)
contour(x,y,FMSPIN)
hold on;
plot(pathap1(:,1), pathap1(:,2),'--bo','Linewidth', 1.5)
plot(pathap2(:,1), pathap2(:,2),'--m+','Linewidth', 1.5)
plot(r1(1),r1(2),'ro','Linewidth',3)
plot(r2(1),r2(2),'ro','Linewidth',3)
plot(r3(1),r3(2),'ro','Linewidth',3)
plot(r4(1),r4(2),'ro','Linewidth',3)
plot(r5(1),r5(2),'ro','Linewidth',3)

xlabel('x')
ylabel('y')
colorbar
title('log(||F_{MSPIN}(x) +1 ||)','Fontsize', 16)
legend('log(||F_{MSPIN}(x) +1 ||)','MSPIN Approx M','MSPIN Approx NM','Real Solution','Location','SouthWest')
hold off;
% det(Japrox(pathap1))

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


function fval=f(b)
x=b(1); y = b(2);
% evaluate function
fval=[((x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -9*y+32)^3- y^3); (x+y + 2)];
end



function fval=f_mspin(b)
x = b(1); y = b(2);
% evaluate function
fval=[(x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -10*y+32); (y^5 - 4*y^4 - 12*y^3 +34*y^2 +11*y -30)];
end


function J_exact = Jex(b)
x = b(1); y = b(2);
J_exact= [3*(x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -9*y+32)^2, (3*(-5*y^4+16*y^3+ 36*y^2 -68 *y -9)*(x-y^5+ 4*y^4 + 12*y^3 -34*y^2 -9*y+32)^2 - 3* y^2); 1 1];
end

function J_mspin = Jmspin(b)
x = b(1); y = b(2);
J_mspin= [1 (-5*y^4+16*y^3+ 36*y^2 -68 *y -10);  0 (5*y^4-16*y^3- 36*y^2 +68 *y +11)] ;
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
