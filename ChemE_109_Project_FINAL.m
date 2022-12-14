%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ChemE 109 Final Project
%         Aaron Lee
%         UID: 505 540 473
%         Steady State Concentration Profile
%         & Time Evolution of Concentration Profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all;
clear cache

tic

% Paramters

toggle1 = true;
toggle2 = true;
plot1 = false;
plot2 = true;

delta = [0, 0.05, 0.05]; % delta
xi = [0, 0.0, 0.03]; % xi
D = 0.1;
beta = 1.5; % beta
gamma = [0.05, 0.02, 0.02]; % gamma
epsilon = [0.0, 0.1, 0.1]; % epsilon
eta = [0.0, 0.05, 0.05]; % eta
theta = [0.0, 0.1, 0.1];

n = 30; % number nodes
dx = 1/n; % x step size
nodes = n*4-2; % number interior nodes

cases = ['A', 'B', 'C'];

%% ------------------ %% Plotting Part 1: Steady State Solutions -----------------------------

if (toggle1)

    for param = 1:3

        if param == 3
            [y, jacob, func] = findSteadyState(n, delta(param), xi(param), D, beta, gamma(param), epsilon(param), eta(param), theta(param));
        else
            y = findSteadyState(n, delta(param), xi(param), D, beta, gamma(param), epsilon(param), eta(param), theta(param));
        end
    
        % Divide vars into 4 species

        ya = y(1:n);
        yb = y(n+1:2*n);
        yu(1) = 0;
        yu(2:n) = y(2*n+1:3*n-1);
        yf(1) = 0;
        yf(2:n) = y(3*n:4*n-2);
        x = linspace(0,1,n);
    
        % Output all nodes at all time, spacesteps

        fprintf('\nSteady State Concentrations: Case(%c)\n',param)
        fprintf('node|       ya|       yb|       yu|       yf|\n')
        for j = 1:1:n
            fprintf('%4i| %7.2e| %7.2e| %7.2e| %7.2e|\n',j, ya(j), yb(j),yu(j),yf(j))
        end
        
    
        %% Plotting Part 1: Steady State Solutions
    
        if (plot1)
            figure(param)
            set(gcf, 'Position', [75 75 1275 750])
            hold on
            grid on
            str = sprintf('Steady State Concentration Case %i', param);
            sgtitle(str,'FontSize',18) 
            %------------------------------------------------
            subplot(2,2,1) 
            plot(x,ya,'color',[1 0 0],'LineWidth',2)
            xlim([0 1])
            ylim([0 max(ya)])
            grid on
            title('[Species A] vs. x','FontSize',14)
            set(gca,'LineWidth',2,'FontSize',12)
            %------------------------------------------------
            subplot(2,2,2) 
            plot(x,yb,'c','LineWidth',2)
            xlim([0 1])
            ylim([0 max(yb)])
            
            grid on
            title('[Species B] vs. x','FontSize',14)
            set(gca,'LineWidth',2,'FontSize',12)
            %------------------------------------------------
            subplot(2,2,3) 
            plot(x,yu,'g','LineWidth',2)
            xlim([0 1])
            if max(yu) == 0
                ylim([0 max(y)])
            else
                ylim([0 max(yu)])
            end
            grid on
            title('[Species U] vs. x','FontSize',14)
            set(gca,'LineWidth',2,'FontSize',12)
            %------------------------------------------------
            subplot(2,2,4) 
            plot(x,yf,'color', [1 0 1],'LineWidth',2)
            xlim([0 1])
            if max(yf) == 0
                ylim([0 max(y)])
            else
                ylim([0 max(yf)])
            end
            grid on
            title('[Species F] vs. x','FontSize',14)
            set(gca,'LineWidth',2,'FontSize',12)  
        end
    end

end


%% -------------------- %% Plotting Part 2: Transient State Solutions --------------------
% using Explicit Euler

if (toggle2)

    dt = 0.001;
    
    % initial conditions
    
    dim = 4*n-2;
    y_old = zeros(dim,1);
    y_new = zeros(dim,1);
    convergence(1) = .01;
    iter = 1;

    INTERVAL = 3; % time interval to plot
    MAX_ITER = 20000;
    limit = 0.99;
    
    figure(4)    
    hold on
    grid on
    xlim([0 1])
    %ylim([0 max(y)])
    set(gca,'LineWidth',2,'FontSize',12)
    set(gcf,'Position',[75 20 750 900])
    xlabel('x')
    ylabel('yi')
    
    %% Plot Part 2

    % EE until converges to steady state
    while (convergence(iter) < limit && iter < MAX_ITER)   
        yy = num2cell(y_old);
        F = func(yy{1:n*4-2}); % gets func from last test case
        
        for i = 1:1:dim    
            y_new(i) = y_old(i) + dt*F(i);
        end
          
        y_old = y_new;
        iter = iter + 1;
        convergence(iter) = min(abs(y_new./y));
        
        time = dt*iter;
        if mod(time, INTERVAL) == 0 % only plot every INTERVAL number of seconds
            
            str = sprintf('Concentration of Species Over Time in %i(s) Increments', INTERVAL);
            title(str,'FontSize',14)
            

            % format data
            yat = y_new(1:n);
            ybt = y_new(n+1:2*n);
            yut(1) = 0;
            yut(2:n) = y_new(2*n+1:3*n-1);
            yft(1) = 0;
     
            yft(2:n) = y_new(3*n:4*n-2);       

            if (plot2)

                grid on
                plot(x,yat,'color', [1 0 0],'LineWidth', 1);
                plot(x,ybt,'c','LineWidth', 1);
                plot(x,yut,'g','LineWidth', 1);
                plot(x,yft,'color', [1 0 1],'LineWidth', 1);
                legend('[Ya]','[Yb]','[Yu]','[Yf]')
    
                %{
                hold on
                grid on
                subplot(2,2,1) 
                title('[Species A] vs. x','FontSize',14)
                plot(x,yat,'color', [1 0 0],'LineWidth', 1);
    
                grid on
                subplot(2,2,2) 
                title('[Species B] vs. x','FontSize',14)
                plot(x,ybt,'c','LineWidth', 1);
    
                grid on
                subplot(2,2,3) 
                title('[Species U] vs. x','FontSize',14)
                plot(x,yut,'g','LineWidth', 1);
    
                grid on
                subplot(2,2,4) 
                title('[Species F] vs. x','FontSize',14)
                plot(x,yft,'color', [1 0 1],'LineWidth', 1);
                %}
            
            end
            
        end
    end


end

toc

%% Function Definitions

function [ynew jacob func] = findSteadyState(n, delta, xi, D, beta, gamma, epsilon, eta,theta)
    % returns steady state values 
    % distribution, general jacobian, and functions for each node
    % n is the number of nodes to use for computation
    % input parameters are coefficients of differential equations
    
    % Parameters
    
    dx = 1/n;
    nodes = n*4-2;
    
    % Create Vectors of Symbolic Variables
    
    ya = sym('ya',[n+2,1]);
    yb = sym('yb',[n+2,1]);
    yu = sym('yu',[n+1,1]);
    yf = sym('yf',[n+1,1]);
    
    %% Boundary Conditions
    
    % @ x = 0
    ya(1) = ya(3) - 2*dx*(ya(2)-1);
    yb(1) = yb(3) - 2*dx*(yb(2)-beta);
    yu(1) = 0;
    yf(1) = 0;
    
    % @ x = 1
    ya(n+2) = ya(n) - 2*dx*epsilon*ya(n+1);
    yb(n+2) = yb(n) - 2*dx*eta*yb(n+1)^2;
    yu(n+1) = yu(n-1) - 2*dx*theta*yu(n);
    yf(n+1) = yf(n-1);
    
    %% Create Symbolic Functions
    % function order is all of ya, then all of yb, etc
    % resulting 4nx4n matrix is ordered with rows 1 to n equating to ya
    % columns 1 to n are ya,i;  col n+1 to 2n are yb,i etc.
    
    for i = 1:nodes
        
        k = i + 1;  
        % k represents node i on x axis
        % % ya(1) is ya0, the computation node      
    
        if i <= n
            f(i) = (D*(ya(k+1) -2*ya(k) + ya(k-1)))/(dx^2) - ya(k)*yb(k)^2 - gamma*ya(k); 
            
        elseif i <= 2*n
            k = k - n; % reorder the index
            f(i) = (D*(yb(k+1) -2*yb(k) + yb(k-1)))/(dx^2) -2*ya(k)*yb(k)^2 - delta*yb(k); 
            
        elseif i <= 3*n-1
            k = k - (2*n);
            f(i) = (D*(yu(k+1) -2*yu(k) + yu(k-1)))/(dx^2) + delta*yb(k+1) - xi*yu(k); 
            
        elseif i <= 4*n-2
            k = k - (3*n-1);
            f(i) = (D*(yf(k+1) -2*yf(k) + yf(k-1)))/(dx^2) + xi*yu(k);    
        end
    end
    
    % Order Node Y Vectors
    
    Y = sym('y', [nodes,1]);
    Y(1:n) =        ya(2:n+1);
    Y(n+1:2*n) =    yb(2:n+1);
    Y(2*n+1:3*n-1) =  yu(2:n);
    Y(3*n:4*n-2) =  yf(2:n);
    
    % Find Jacobian
    
    jacob = jacobian(f,Y);

    %% Newtons method for Jdy = F
    % initial guess
    y = zeros(nodes,1);
    
    % converting symbolic functions into matlab functions faster runtime

    jacob = matlabFunction(jacob, 'vars', Y);
    func = matlabFunction(f, 'vars', Y);
    iter = 1;
    error(1) = 100;
    ERROR = 1e-5;
    MAX_ITER = 50;

    % Multivar Newton's method
    while iter < MAX_ITER && error(iter) > ERROR

        yold = num2cell(y);
        J = jacob(yold{1:nodes});
        F = func(yold{1:nodes});
        F = reshape(F, [nodes 1]);
        ynew = y - J\F;
        y = ynew;
        iter = iter + 1;
        error(iter) = max(norm(F));
    end

end


