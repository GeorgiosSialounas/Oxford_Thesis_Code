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
npts=201; %npts=2001 gives good resolution
X=linspace(-1,3,npts);
Y=linspace(-1,3,npts);
r1=[1  ;1];
set(gca,'Fontsize',14)
tol=1e-10;



[x,y]=meshgrid(X,X); % set up a grid of starting guesses for Newton's method
F = zeros(size(x));
FMSPIN = zeros(size(x));
for  i = 1 : npts
    for j = 1 : npts
        F(i,j) = log(norm(f([x(i,j),y(i,j)]),2)+1);
        FMSPIN(i,j) = log(norm(f_mspin([x(i,j),y(i,j)]),2)+1);
   
    end
end
[~,~,pathh,~]=mynewt_mspinapprox_path(-0.86,2.85, tol);
[~,~,pathhgood,~]=mynewt_mspinapproxgood_path(-0.86,2.85, tol);
[~,~,pathnewt]=mynewt_path(-0.86,2.85, tol);
figure(6);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
subplot(1,2,1)
hold on;
contour(x,y,F)
% plot(pathnewt(:,1), pathnewt(:,2),'--b*','Linewidth', 1.5)
hold off;

xlabel('x')
ylabel('y')
colorbar
title('log(||F(x) +1 ||)')
subplot(1,2,2)
contour(x,y,FMSPIN)
hold on;
plot(pathhgood(:,1), pathhgood(:,2),'--m+','Linewidth', 1.5)
plot(pathh(:,1), pathh(:,2),'--bo','Linewidth', 1.5)
plot(1,1,'ro')
xlabel('x')
ylabel('y')
colorbar
title('log(||F_{MSPIN}(x) +1 ||)')
legend('log(||F_{MSPIN}(x) +1 ||)','MSPIN Approx. M','MSPIN Approx. NM','Real Solution','Location','SouthEast')
hold off;
% define the three roots of f(z)

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
cp=zeros(size(x));
cit=zeros(size(x));
sing_mat = zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k,sing]=mynewt(x(i,j),y(i,j),tol);
        %         norm(res-r1)
        cit(i,j) = k;
        sing_mat(i,j) = sing;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
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

% subplot(3,3,7)
% surf(x,y,sing_mat), view(2), shading interp, axis equal tight
% title('Original Jacobian')


% Newton with exact MSPIN jacobian
cp=zeros(size(x));
cit=zeros(size(x));
sing_mat = zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res, k, sing]=mynewt_mspinexact(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        sing_mat(i,j) = sing;
        %         %         norm(res-r1)
        %         r1
        %         res
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
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

% subplot(3,3,8)
% surf(x,y,sing_mat), view(2), shading interp, axis equal tight
% title('True MSPIN Jacobian')


% Newton with approximate MSPIN jacobian
cp=zeros(size(x));
cit=zeros(size(x));
sing_mat = zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        
        [res, k,sing]=mynewt_mspinaprox(x(i,j),y(i,j),tol);
        cit(i,j) = k;
        sing_mat(i,j) = sing;
        if norm(res-r1)<tol       % if there is convergence to root 1
            cp(i,j)=1;            % then record this
            cmm(i,j,:)=[1,0,0];   % record colour to plot as red
            
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
% subplot(3,3,9)
% surf(x,y,sing_mat), view(2), shading interp, axis equal tight
% title('Approx. MSPIN Jacobian')




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
% subplot(3,3,9)
% surf(x,y,sing_mat), view(2), shading interp, axis equal tight
% title('Approx. MSPIN Jacobian')

[res, k,path]=mynewt_mspinexact_path(-0.86,0,tol);
hold off
figure;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
plot(path(:,1),path(:,2),'-b*')
xlabel('x')
ylabel('y')
title('path plot')

[~, ~,path,~]=mynewt_mspinapprox_path(-0.86,2.85, tol);
[~, ~,pathgood,~]=mynewt_mspinapproxgood_path(-0.86,2.85, tol);
hold off
figure;

xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
subplot(1,2,1)
plot(pathgood(:,1),pathgood(:,2),'-r+')

hold on;

plot(path(:,1),path(:,2),'-b*')
legend('MSPIN Approx. M','MSPIN Approx. NM','SouthWest')

plot(1,1,'ko','Linewidth',1)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
xlabel('x')
ylabel('y')
title('Approximate MSPIN Jacobians: sequence of iterates')
grid on
hold off

dets  = zeros(1, length(path));
condnums  = zeros(1, length(path));

eigvals = zeros(size(path));

for i = 1 : length(path)
    Jma = Japrox2(path(i,:));
    Jma = Japrox2(path(i,:));
    Jmainv = Jma\eye(2);
    dets(i) = det(Jma);
    condnums(i) = cond(Jma);
    eigvals(i,:) = eig(Jma);
    mineigvals(i,:) = min(eig(Jma));
    detsinv(i) = det(Jmainv);
    eigvalsinv(i,:) = eig(Jmainv);
    
end

% figure;
subplot(1,2,2)
plot(1:length(condnums), condnums, '--b*')% plot(1:length(detsinv), detsinv, '--b*')
xlabel('Iterate')
ylabel('Condition Number')
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
title('MSPIN Approx. NM: Cond. Number')
grid on
hold off;
function [xnew,k, sing]=mynewt(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f(xold);
k=0;
sing = 0;
while (norm(fxold,Inf)>tol)
    Je = Jex(xold);
    
    if det(Je)==0
        sing = 1;
    end
    xnew=xold-(Je)\fxold;
    xold=xnew;
    fxold=f(xold);
    k=k+1;
    if k >50
        break;
    end
end
xnew=xold;
end



% function for original jacobian
function [xnew,k,path]=mynewt_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f(xold);
k=0;
path = [];
while (norm(fxold,Inf)>tol)
    
    path(k+1,1:2)= xold';
    xnew=xold-(Jmspin(xold))\fxold;
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

% function for good approximate mspin jacobian
function [xnew,k,path, sing]=mynewt_mspinapprox_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
path = [];
sing = 0;
while (norm(fxold,Inf)>tol)
    if det(Japrox2(xold))==0
        sing = 1;
    end
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

% function for good approximate mspin jacobian
function [xnew,k,path, sing]=mynewt_mspinapproxgood_path(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
path = [];
sing = 0;
while (norm(fxold,Inf)>tol)
    if det(Japrox(xold))==0
        sing = 1;
    end
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



% function for bad approximate jacobian
function [xnew,k,sing]=mynewt_mspinexact(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
sing = 0 ;
while (norm(fxold,Inf)>tol)
    Jm = Jmspin(xold);
    if det(Jm)==0
        sing = 1;
    end
    xnew=xold-(Jm)\fxold;
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
function [xnew,k,sing]=mynewt_mspinaprox(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
sing = 0;
while (norm(fxold,Inf)>tol)
    Jap = Japrox(xold);
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

function fval=f(x)
% evaluate function
fval=[((x(1)-x(2)^3 +1)^3 -x(2)^3); (2*x(1)+3*x(2)-5)];
end

function fval=f_mspin(x)
% evaluate function
fval=[((x(1)-x(2)^3 +1) -x(2));( (2/3)*x(2)^3 + (5/3)*x(2)-7/3)];
end


function J_exact = Jex(x)
J_exact= [3*(x(1)-x(2)^3+1)^2 , ((3*(x(1)-x(2)^3+1)^2)*(-3*x(2)^2)-3*x(2)^2); 2 3];
end

function J_mspin = Jmspin(x)
J_mspin= [1 (-3*x(2)^2 -1);  0 (2*x(2)^2+5/3)] ;
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
