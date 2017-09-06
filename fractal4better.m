clc;
clear;
close all;

% G=F1(x1,x2)=(x1-x2^3-1)^3-x2^3
% H=F2(x1,x2)=x1+2*0.5*exp(x2)-1-2*0.5
% The exact solution is x_exact=[1,1]';

% g=x1-x2^3-x2-1
% h= x2- log(1-(1/(2*mu))*(x2^3+x2)) %0.5*x2^3+1.5*x2-2

mu=0.5;
tol=1e-9;
npts=201; %npts=2001 gives good resolution
X=linspace(-1,3,npts);
Y = linspace(-2,2,npts);
[x,y]=meshgrid(X,Y); % set up a grid of starting guesses for Newton's method

% define the two roots of f(z)
r1=[1; 0];

cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt(x(i,j),y(i,j),tol);
%         res
        %         norm(res-r1)
         cit(i,j) = k;      
%         cit(i,j) = 50;
        if norm(res-r1)<tol       % if there is convergence to root 1
             
%             cit(i,j) = k;
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

figure(1);

hold on;
subplot(2,2,1)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
% colorbar
title('Original Jacobian','Fontsize', 14)
cpnewt = cp;
figure(2)

hold on;
subplot(2,2,1)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Original Jacobian','Fontsize', 14)
colorbar


figure;
subplot(1,2,1)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
surf(x,y,cp,cmm), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
% colorbar
title('Original Jacobian','Fontsize', 14)
cpnewt = cp;
subplot(1,2,2)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
surf(x,y,cit), view(2), shading interp, axis equal tight
xlabel('x')
ylabel('y')
title('Original Jacobian','Fontsize', 14)
colorbar


% Exact MSPIN
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt_mspinexact(x(i,j),y(i,j),tol);
        %         cit(i,j) = k;
%         cit(i,j) = k;
cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
%             cit(i,j) = k;
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

% Good Aprox MSPIN
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
        %         cit(i,j) = k;
%         cit(i,j) = k;
cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
%             cit(i,j) = k;
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



% Bad Aprox MSPIN
cp=zeros(size(x));
cit=zeros(size(x));
% run Newton's method for each starting guess
for i=1:npts
    %i
    for j=1:npts
        [res,k]=mynewt_mspinaprox2(x(i,j),y(i,j),tol);
        %         norm(res-r1)
        %         disp('here')
        %         i
        %         j
        %         cit(i,j) = k;
%         cit(i,j) = k;
 cit(i,j) = k;
        if norm(res-r1)<tol       % if there is convergence to root 1
%             cit(i,j) = k;
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



[res_mspinex, k,path_mspinex]=mynewt_mspinexact_path(-1.5,1.5,tol);
[res_newt, k,path_newt]=mynewt_path(-1.5,1.5,tol);
[res_mspinap1, k,path_mspinap1]=mynewt_mspinaprox1_path(-1.5,1.5,tol);
[res_mspinap2, k,path_mspinap2]=mynewt_mspinaprox2_path(-1.5,1.5,tol);

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



% newton path
[~,~,path_newt1]=mynewt_path(1.5,-1.5,tol);
[~,~,path_newt2]=mynewt_path(1.5,1.5,tol);
[~,~,path_newt3]=mynewt_path(-1.5,-1.5,tol);
[~,~,path_newt4]=mynewt_path(-1.5,1.5,tol);
% MSPIN_exact path
[~,~,path_ex1]=mynewt_mspinexact_path(1.5,-1.5,tol);
[~,~,path_ex2]=mynewt_mspinexact_path(1.5,1.5,tol);
[~,~,path_ex3]=mynewt_mspinexact_path(-1.5,-1.5,tol);
[~,~,path_ex4]=mynewt_mspinexact_path(-1.5,1.5,tol);
% MSPIN approx1 path
[~,~,path_ap11]=mynewt_mspinaprox1_path(1.5,-1.5,tol);
[~,~,path_ap12]=mynewt_mspinaprox1_path(1.5,1.5,tol);
[~,~,path_ap13]=mynewt_mspinaprox1_path(-1.5,-1.5,tol);
[~,~,path_ap14]=mynewt_mspinaprox1_path(-1.5,1.5,tol);
% MSPIN approx2 path
[~,~,path_ap21]=mynewt_mspinaprox2_path(1.5,-1.5,tol);
[~,~,path_ap22]=mynewt_mspinaprox2_path(1.5,1.5,tol);
[~,~,path_ap23]=mynewt_mspinaprox2_path(-1.5,-1.5,tol);
[~,~,path_ap24]=mynewt_mspinaprox2_path(-1.5,1.5,tol);
figure;
hold on;
plot(path_newt1(:,1),path_newt1(:,2),'b--*');
plot(path_newt2(:,1),path_newt2(:,2),'k--+');
plot(path_newt3(:,1),path_newt3(:,2),'m--o');
plot(path_newt4(:,1),path_newt4(:,2),'g--+');
plot(r1(1),r1(2),'ro','LineWidth',1.5)
hold off;
legend('(1.5,-1.5)','(1.5,1.5)','(-1.5,-1.5)','(1.5,-1.5)','Solution')
title('Paths using Newton Method','Fontsize',16)
xlabel('x')
ylabel('y')


path1 = path_newt1;
path2 = path_newt2;
path3 = path_newt3;
path4 = path_newt4;
dets  = zeros(1, length(path));
% condnums  = zeros(1, length(path));

% eigvals = zeros(length(path));
%%
for i = 1 : length(path1)
    
    Jma = Jex(path1(i,:));
    dets(i) = det(Jma);
    condnums(i) = cond(Jma);
    eigvals(i,:) = eig(Jma)';
    mineigvals(i,:) = min(eig(Jma));
end
for i = 1 : length(path2)
        Jma2 = Jex(path2(i,:));
    dets2(i) = det(Jma2);
    condnums2(i) = cond(Jma2);
    eigvals2(i,:) = eig(Jma2)';
    mineigvals2(i,:) = min(eig(Jma2));
end
for i = 1 : length(path3)
        Jma3 = Jex(path3(i,:));
    dets3(i) = det(Jma3);
    condnums3(i) = cond(Jma3);
    eigvals3(i,:) = eig(Jma3)';
    mineigvals3(i,:) = min(eig(Jma3));
end
for i = 1 : length(path4)
        Jma4 = Jex(path4(i,:));
    dets4(i) = det(Jma4);
    condnums4(i) = cond(Jma4);
    eigvals4(i,:) = eig(Jma4)';
    mineigvals4(i,:) = min(eig(Jma4));
    
end
%%
figure;
% subplot(1,2,2)
plot(1:length(condnums), condnums, '--b*')% plot(1:length(detsinv), detsinv, '--b*')
hold on;
plot(1:length(condnums2), condnums2, '--r*')% plot(1:length(detsinv), detsinv, '--b*')
plot(1:length(condnums3), condnums3, '--m*')% plot(1:length(detsinv), detsinv, '--b*')
plot(1:length(condnums4), condnums4, '--k*')% plot(1:length(detsinv), detsinv, '--b*')
xlabel('Iterate')
ylabel('Condition Number')
legend('(1.5,-1.5)','(1.5,1.5)','(-1.5,-1.5)','(-1.5,1.5)','Location','NorthWest')
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
title('Newton Method: Condition Number','Fontsize',16)
grid on

path1 = path_newt1;
path2 = path_newt2;
path3 = path_newt3;
path4 = path_newt4;
dets  = zeros(1, length(path));
% condnums  = zeros(1, length(path));



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


% function for bad mspin approx
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




function [xnew,k]=mynewt(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f(xold);
k=0;
while (norm(xold-[1,0]',Inf)>tol)
    xnew=xold-(Jex(xold))\fxold;
    xold=xnew;
    fxold=f(xold);
    k=k+1;
    if k >1000
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

% function for approximate mspin jacobian
function [xnew,k]=mynewt_mspinaprox2(x,y,tol)
% main Newton loop
xold=[x,y]';
fxold=f_mspin(xold);
k=0;
while (norm(fxold,Inf)>tol)
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

function fval=f(x)
% evaluate function
mu=0.5;
fval=[((x(1)-x(2)^3-1)^3-x(2)^3); (x(1)+2*mu*exp(x(2))-1-2*mu)];
end

function fval=f_mspin(x)
mu=0.5;
% evaluate function
fval=[(x(1)-x(2)-x(2)^3-1); (x(2)- log(1-1/(2*mu)*(x(2)^3 + x(2)) ))];
end



function J_exact = Jex(x)
mu=0.5;
J_exact= [3*(x(1)-x(2)^3-1)^2, (-9*x(2)^2*(x(1)-x(2)^3-1)^2 - 3*x(2)^2); 1, 2*mu*exp(x(2))];
end

function J_mspin = Jmspin(b)
mu=0.5;
x = b(1); y = b(2);
J_mspin= [1 (-3*y^2-1); 0 ( 1+(0.5/mu)*(3*y^2 + 1)/(1-(1/(2*mu))*(y^3 + y) ) )];
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

function J_mspin_aprox = Japrox2(x)

J = Jex(x);
L = tril(J);
J_mspin_aprox = L\J;
end



