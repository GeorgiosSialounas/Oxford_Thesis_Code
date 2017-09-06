function MSPIN_codes
%% This is a simple example with two unknowns

% G=F1(x1,x2)=(x1-x2^3+1)^3-x2^3
% H=F2(x1,x2)=x1+2x2-3
% The exact solution is x_exact=[1,1]';

% g=x1-x2^3-x2+1
% h=0.5*x2^3+1.5*x2-2

% Written by Lulu Liu
% Extreme Computing Research Center
% King Abdullah University of Science and Technology (KAUST)
% Nov, 2014

%%
clc;clf;
format long;

%% the exact solution
x_exact=[1,1]';

%% Initialization
x0=[0.0,0.0]';      % initial guess

max_its=50;         % maximum iterations for nonlinear solvers
tol=1.e-6;          % tolerance for nonlinear solvers
ksp_tol=1.e-8;      % tolerance for linear solvers

% options
linesearch = true;  % using linesearch
ls_monitor=true;    % linesearch monitor
snes_monitor=true;  % nonlinear solver monitor
plot_path = true;   % plot paths of approximated solutions

%%
fprintf('======================================================\n');
fprintf('                      MSPIN                           \n');
fprintf('======================================================\n');
x_soln=x0;
F=fun_MSPIN(x0);
normF0=norm(F);
iter=0;
stp=1;

%Records for solutions and residuals
XX_MSPIN=[]; % solutions 
R_MSPIN=[];  % Residuals of the preconditioned system
RR_MSPIN=[]; % Residuals of the original system

XX_MSPIN(1,:)=x0';
R_MSPIN(1,:)=F';
FX=fun_old(x_soln);
RR_MSPIN(1,:)=FX';

while norm(F)/normF0>tol
    iter=iter+1;
    
   %% Solve for Newton direction d using analytical Jacobian(Matrix-vector)
   % I implemented (2.28) in my paper, but you can comment the following
   % lines and try the next block (approximate Jacobian)
   
   [d,flg]=gmres(@(x)Jv_MSPIN(x,x_soln),-F,[],ksp_tol);
   if(flg>0)
       error('The linear solver does not converge at step %d\n',iter);
   end 

   %% Solve for Newton direction d using approximate Jacobian
   %  Adding the following lines for you 
   %  You will find the Jacobian of the original system will be singular at
   %  pv=(p,v), but MSPIN works well at x_soln=(u,v)
   
    % caculate g and h gh=[g;h]
      gh=fun_MSPIN(x_soln);
    % caculate [p;v]=[u-g;v] (x_soln=[u,v])
      pv=x_soln-[gh(1);0];    
      
    % Jacobian is singular at [p;v] starting from(0,0) (the first row is zeros), 
    % that is why I choose analytical Jacobian in my paper
      %J=Jac_old(pv) 
      
    % Caculate approximate Jacobian at (u,v) instead of (p,v) for the formula (2.29) in my paper  
      J=Jac_old(x_soln); 
      L=tril(J);    
      J=L\J; % L^{-1}J
      d=-J\F;

   %% Linesearch        
      if(linesearch)
       % fprintf('Iter %d, Linesearch begins...........\n',iter);
        F=fun_MSPIN(x_soln);        
        f=0.5*F'*F;
        dg=-2*f;
        stp=1;       
        stp = ls( x_soln,d, dg, f, stp, [],@(x)fun_MSPIN(x)); 
        if(ls_monitor)
            fprintf('Iter %d   Line search: stepsize = %12.9e \n',iter-1,stp);
        end
      end    
    
   %% update solution   
    x_soln=x_soln+stp*d;
    XX_MSPIN(iter+1,:)=x_soln';
    F=fun_MSPIN(x_soln);
    R_MSPIN(iter+1,:)=F';
    FX=fun_old(x_soln);
    RR_MSPIN(iter+1,:)=FX';    
end

if(snes_monitor)
    fprintf('# MSPIN: the Jacobian is analytical form here.\n');
    fprintf('The intial value is x = (%d, %d)\n',x0(1),x0(2));
    fprintf('-------------------------------------------------------------------------------\n');     
    fprintf('step      x=(x1,x2)              ||F(x)||          ||T(x)||         T1(x1,x2)       T2(x1,x2)\n');
    fprintf('----   ------------------     --------------   --------------   ------------    ------------\n');
    
    for i=1:iter+1    
        fprintf('%3i   [%f, %f]      %e      %e      %e     %e\n',i-1,XX_MSPIN(i,1),...
            XX_MSPIN(i,2),norm(RR_MSPIN(i,:)),norm(R_MSPIN(i,:)),R_MSPIN(i,1),R_MSPIN(i,2));
    end
end
fprintf('||x-x*|| = %e   ||F(x)|| = %e\n',norm(x_soln-x_exact),norm(R_MSPIN(end,:)));
x_soln


figure(1)
for i=1:iter+1
normF_MSPIN(i)=norm(R_MSPIN(i,:));
end
semilogy(1:iter+1,normF_MSPIN,'m-o');
title('Residual history of MSPIN')

end
%% Function Evaluation
function F=fun_old(x)
F=zeros(2,1);
F(1)=(x(1)-x(2)^3+1)^3-x(2)^3;
F(2)=x(1)+2*x(2)-3;
end



function F=fun_MSPIN(x)
F=zeros(2,1);
F(1)=x(1)-x(2)-x(2)^3+1;
F(2)=0.5*x(2)^3+1.5*x(2)-2;
end
%% Jac-vec multiplication
function y=Jv_MSPIN(x,r)
y=ones(2,1);
y(1)=x(1)-(1+3*r(2)*r(2))*x(2);
y(2)=1.5*(1+r(2)*r(2))*x(2);
end


%% Jacobian
function J=Jac_old(x)
J=[3*(x(1)-x(2)^3+1)^2, -9*x(2)^2*(x(1)-x(2)^3+1)^2-3*x(2)^2;1,2];
end


