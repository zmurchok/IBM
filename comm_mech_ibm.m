%Communication Mechanisms Enrich Pattern Formation in an Individual-Based Model of Collective Behaviour
%Cole Zmurchok and Gerda de Vries
%2017

%This .m file contains an implementation of the IBM described in the paper. Running the file as is will produce figures similar to those in Figure 4 and 5 of the paper.

%The file contains the following functions:
%main(...) which calls the IBM and the plotting functions.
%ibm(...) which implements the IBM as described in the paper.
%plotter(...) which produces plots and implements the kernel smoothing density estimation
%ibm_density_dependent_speed(...) which implements the IBM with density-dependent movement speed

N = 50; %number of individuals
steps = 3000; %number of timesteps
plotmin = 2000; %plotmin < t < plotmax gives the values of time that will be plotted
plotmax = 3000; %the default setting is plot only the last 1000 timesteps

%stat1 will produce a stationary pulse and save a file called stat1.mat with the workspace variables for reference
%main() requires the followign parameters as input
%title and filename: '(a) Stationary pulse 1'
%submodel 'M1'
%the remaining parameters are: lambda_1, lambda_2, qr, qal, qa, N, steps, plotmin, plotmax

stat1 = main('(a) Stationary pulse 1','M1',0.2,0.9,2.4,0,2,N,steps,plotmin,plotmax);
clear stat1;

stat2 = main('(b) Stationary pulse 2','M2',0.2,0.9,0.5,0,4,N,steps,plotmin,plotmax);
clear stat2;

ripples = main('(c) Ripples','M5',0.2,0.9,1.1,2,1.5,N,steps,plotmin,plotmax);
clear ripples;

feathers = main('(d) Feathers','M3',0.2,0.9,6.4,0,6,N,steps,plotmin, ...
               plotmax);
clear feathers;

pulse = main('(e) Traveling pulse','M1',0.2,0.9,0.5,2,1.6,N,steps,plotmin,plotmax);
clear pulse;

train = main('(f) Train','M3',6.67,30,0,2,0,N,steps,plotmin,plotmax);
clear train;

zigzag = main('(g) Zigzag pulse','M2',0.2,0.9,1,2,6,N,steps,plotmin,plotmax);
clear zigzag;

breather = main('(h) Breathers','M4',0.2,0.9,1,0,2,N,steps,plotmin,plotmax);
clear breather;

travbreather = main('(i) Traveling breather','M4',0.2,0.9,4,2,4,N,steps,plotmin, ...
                   plotmax);
clear travbreather;

function y = main(filename,submodel,lambda1,lambda2,qr,qal,qa,N,steps,plotmin,plotmax)
    y = ibm(submodel,lambda1,lambda2,qr,qal,qa,N,steps);
    plotter(filename,y,N,plotmin,plotmax);
    save(filename,'y');
end

function y = plotter(filename,input,N,a,b)
    %this function produces *nice* plots! and implements the kernel estimate of the density
    %See for example: https://en.wikipedia.org/wiki/Kernel_density_estimation
    fact = 1;
    width=2.5*fact;
    height=2.5*fact;
    x0 = 5;
    y0 = 5;
    fontsize = 10;

    filename1 = strcat(filename,'.eps');
    filename2 = strcat(filename,'.tif');
    filename3 = strcat(filename,'_density.eps');
    filename4 = strcat(filename,'_density.tif');
    filename5 = strcat(filename,'_density_colorbar.eps');
    filename6 = strcat(filename,'_density_colorbar.tif');

    domainsize = 10;
    x = a:b;
    vector = input;
    figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto')
    hold on;
    xlabel({'$x$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    ylabel({'$t$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    title(filename,'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    for i=1:N
        scatter(vector(x,i),x,0.1,'.');
    end
    xlim([0 10]);
    print(1,filename1,'-depsc','-r300')
    print(1,filename2,'-dtiff','-r300')
    close all
    %density estimation:
    gx = (domainsize/N):(domainsize/N):domainsize; %grid for estimation
    phi=@(x)(exp(-.5*x.^2)/sqrt(2*pi)); %standard normal
    h = 0.1;
    for t=a:b
        %kernel density estimation loop
        %difference is a fancy way of calculating x_i-gx_i for all i
        difference = diag(vector(t,:))*ones(N,N) - ones(N,N)*diag(gx);
        d(t,:) = sum(phi(difference/h))/N/h;
    end
    figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto')
    hold on;
    xlabel({'$x$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    ylabel({'$t$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    title(strcat(filename,' Density'),'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    %make top-down surface plot:
    mesh(gx,x,d(x,:));
    view(2)

    print(1,filename3,'-depsc','-r300')
    print(1,filename4,'-dtiff','-r300')
    close all

    figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto')
    hold on;
    xlabel({'$x$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    ylabel({'$t$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    title(strcat(filename,' Density'),'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',fontsize,'FontName','Times')
    mesh(gx,x,d(x,:));
    view(2)
    colorbar;

    print(1,filename5,'-depsc','-r300')
    print(1,filename6,'-dtiff','-r300')
    close all
end


function y=ibm(submodel,lambda1,lambda2,qr,qal,qa,N,steps)

    %Model parameters (fixed)
    L = 10;
    y0 = 2;
    sr = 0.25;
    sal = 0.5;
    sa = 1;
    mr = sr/8;
    mal = sal/8;
    ma = sa/8;

    %Numerical parameters (fixed)
    gamma = 0.1;
    dt = 0.05;
    A = 2;
    normalization = A*L/N;

    %The following if statement deals with the varying
    %submodels. A vector of flags is created, which corresponds to
    %which groups of individuals are included in the nonlocal
    %social interaction integral terms from the Eulerian model. The
    %flags vector follows the following format:
    %flags=[rmr_r,lmr_r,rml_r,lml_r,rmr_al,lmr_al,rml_al,lml_al,rmr_a,lmr_a,rml_a,lml_a],
    %where rmr is right-moving right-individuals i.e. u^+(x+s)
    %and lmr is left-moving right-individuals i.e. u^-(x+s). The r,
    %al, and a represent repulsion, aligment and attraction,
    %respectively.

    %The nonlocal interaction terms in submodels 1,2,4 are
    %symmetric for right and left moving individuals, hence we need
    %only to calculate the flag vector for right_moving individuals
    %and multiply by -1 for the left_moving individuals. It should
    %be noted that the q_i (strengths of social interactions) are
    %already incorporated into these flag vectors.

    %We can read the nonlocal social interaction terms from the
    %flag vectors as follows. Look at submodel M1. We have
    %[1,1,-1,-1,0,1,-1,0,1,1,-1,-1]. This means that in the
    %repulsion interaction term, we add rmr + lmr individuals and
    %subtract rml and subtract lml individuals. Indeed, this agrees
    %with the integral from the original model y_r^+ =
    %\int_0^\infty K_r(s) (u(x+s)-u(x-s)) ds, where u = u^+ + u^-.

    if strcmp(submodel,'M1')
        right_moving_flags = [1,1,-1,-1,0,1,-1,0,1,1,-1,-1];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M2')
        right_moving_flags = [1,1,-1,-1,-1,1,-1,1,1,1,-1,-1];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M3')
        right_moving_flags = [1,1,0,0,-1,1,0,0,1,1,0,0];
        left_moving_flags = [0,0,1,1,0,0,1,-1,0,0,1,1];
    elseif strcmp(submodel,'M4')
        right_moving_flags = [0,1,-1,0,0,1,-1,0,0,1,-1,0];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M5')
        right_moving_flags = [0,1,0,0,0,1,0,0,0,1,0,0];
        left_moving_flags = [0,0,1,0,0,0,1,0,0,0,1,0];
    end

    %Randomly assign positions between 0 and domainsize and velocities of
    %+-1
    x = rand(1,N)*L;
    v = 1 - 2.*floor(2.*rand(1,N));

    %indicator function
    chi = @(s) (s>0);
    gaussian = @(a,b,c) exp((-(a-b).^2)/(2*c^2))/sqrt(2*pi*c^2);
    f = @(d) (0.5+0.5*tanh(d-y0));
    %intital storage of the trajectories
    traj = x;

    %Enter time step loop
    for k = 1:steps
        %calculate matrix x_i - x_j. i = row, j = col.
        diffx = diag(x)*ones(N,N);
        diffx = diffx - diffx';
        %take care of periodic boundary conditions, and insert 1000000 into the
        %diagonals so that the individual at x_i doesn't influence itself
        diffxp = mod(diffx,L) + 1e6*eye(N);
        diffxm = -mod(diffx,-L) + 1e6*eye(N);

        %right movers and left movers
        %stands for velocity plus or velocity minus
        vp = (v>0);
        vm = (v<0);

        %calculation of interaction terms for all individuals and all directions
        rmr_r = qr*vp*(gaussian(diffxp,sr,mr).*chi(2*sr-diffxp));
        lmr_r = qr*vm*(gaussian(diffxp,sr,mr).*chi(2*sr-diffxp));
        rml_r = qr*vp*(gaussian(diffxm,sr,mr).*chi(2*sr-diffxm));
        lml_r = qr*vm*(gaussian(diffxm,sr,mr).*chi(2*sr-diffxm));
        %For the revised repulsion kernel:
        % rmr_r = qr*vp*(gaussian(diffxp,0,4*mr).*chi(2*sr-diffxp));
        % lmr_r = qr*vm*(gaussian(diffxp,0,4*mr).*chi(2*sr-diffxp));
        % rml_r = qr*vp*(gaussian(diffxm,0,4*mr).*chi(2*sr-diffxm));
        % lml_r = qr*vm*(gaussian(diffxm,0,4*mr).*chi(2*sr-diffxm));

        rmr_al = qal*vp*(gaussian(diffxp,sal,mal).*chi(2*sal-diffxp));
        lmr_al = qal*vm*(gaussian(diffxp,sal,mal).*chi(2*sal-diffxp));
        rml_al = qal*vp*(gaussian(diffxm,sal,mal).*chi(2*sal-diffxm));
        lml_al = qal*vm*(gaussian(diffxm,sal,mal).*chi(2*sal-diffxm));
        rmr_a = -qa.*vp*(gaussian(diffxp,sa,ma).*chi(2*sa-diffxp));
        lmr_a = -qa.*vm*(gaussian(diffxp,sa,ma).*chi(2*sa-diffxp));
        rml_a = -qa.*vp*(gaussian(diffxm,sa,ma).*chi(2*sa-diffxm));
        lml_a = -qa.*vm*(gaussian(diffxm,sa,ma).*chi(2*sa-diffxm));

        %create output matrix;
        output = vertcat(rmr_r,lmr_r,rml_r,lml_r,rmr_al,lmr_al, ...
                         rml_al,lml_al,rmr_a,lmr_a,rml_a,lml_a);
        %apply submodel by using flags
        signal_p = right_moving_flags*output;
        signal_m = left_moving_flags*output;

        %Calculate the turning rates
        lambda_p = lambda1 + lambda2*f(normalization.*signal_p);
        lambda_m = lambda1 + lambda2*f(normalization.*signal_m);

        %update positions
        x = x + gamma*dt*v;
        %Fix up periodic boundaries in case individuals moved too far
        x = mod(x,L);
        %store positions to plot later
        traj = [traj; x];

        %update velocities as per Eq 4
        rp = rand(1,N);
        rm = rand(1,N);
        v = chi(v).*sign(rp-lambda_p*dt)-chi(-v).*sign(rm-lambda_m*dt);
    end
    y = traj;
end

function y=ibm_density_dependent_speed(submodel,lambda1,lambda2,qr,qal,qa,N,steps)

    %Model parameters (fixed)
    L = 10;
    y0 = 2;
    sr = 0.25;
    sal = 0.5;
    sa = 1;
    mr = sr/8;
    mal = sal/8;
    ma = sa/8;

    %Numerical parameters (fixed)
    gamma = 0.1;
    dt = 0.05;
    A = 2;
    normalization = A*L/N;

    %The following if statement deals with the varying
    %submodels. A vector of flags is created, which corresponds to
    %which groups of individuals are included in the nonlocal
    %social interaction integral terms from the Eulerian model. The
    %flags vector follows the following format:
    %flags=[rmr_r,lmr_r,rml_r,lml_r,rmr_al,lmr_al,rml_al,lml_al,rmr_a,lmr_a,rml_a,lml_a],
    %where rmr is right-moving right-individuals i.e. u^+(x+s)
    %and lmr is left-moving right-individuals i.e. u^-(x+s). The r,
    %al, and a represent repulsion, aligment and attraction,
    %respectively.

    %The nonlocal interaction terms in submodels 1,2,4 are
    %symmetric for right and left moving individuals, hence we need
    %only to calculate the flag vector for right_moving individuals
    %and multiply by -1 for the left_moving individuals. It should
    %be noted that the q_i (strengths of social interactions) are
    %already incorporated into these flag vectors.

    %We can read the nonlocal social interaction terms from the
    %flag vectors as follows. Look at submodel M1. We have
    %[1,1,-1,-1,0,1,-1,0,1,1,-1,-1]. This means that in the
    %repulsion interaction term, we add rmr + lmr individuals and
    %subtract rml and subtract lml individuals. Indeed, this agrees
    %with the integral from the original model y_r^+ =
    %\int_0^\infty K_r(s) (u(x+s)-u(x-s)) ds, where u = u^+ + u^-.

    if strcmp(submodel,'M1')
        right_moving_flags = [1,1,-1,-1,0,1,-1,0,1,1,-1,-1];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M2')
        right_moving_flags = [1,1,-1,-1,-1,1,-1,1,1,1,-1,-1];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M3')
        right_moving_flags = [1,1,0,0,-1,1,0,0,1,1,0,0];
        left_moving_flags = [0,0,1,1,0,0,1,-1,0,0,1,1];
    elseif strcmp(submodel,'M4')
        right_moving_flags = [0,1,-1,0,0,1,-1,0,0,1,-1,0];
        left_moving_flags = -1.*right_moving_flags;
    elseif strcmp(submodel,'M5')
        right_moving_flags = [0,1,0,0,0,1,0,0,0,1,0,0];
        left_moving_flags = [0,0,1,0,0,0,1,0,0,0,1,0];
    end


    %Randomly assign positions between 0 and domainsize and velocities of
    %+-1
    x = rand(1,N)*L;
    v = 1 - 2.*floor(2.*rand(1,N));

    %indicator function
    chi = @(s) (s>0);
    gaussian = @(a,b,c) exp((-(a-b).^2)/(2*c^2))/sqrt(2*pi*c^2);
    f = @(d) (0.5+0.5*tanh(d-y0));
    g = @(s) 1 + tanh(s);

    %storage of the trajectories
    traj = x;

    %Enter time step loop
    for k = 1:steps

        %%%%%%%%%%%%%%%%%%%%%%% TURNING half-STEP

        %calculate matrix x_i - x_j. i = row, j = col.
        diffx = diag(x)*ones(N,N);
        diffx = diffx - diffx';
        diffxp = mod(diffx,L) + 1e6*eye(N);
        diffxm = -mod(diffx,-L) + 1e6*eye(N);

        %right movers and left movers
        %stands for velocity plus or velocity minus
        vp = (v>0);
        vm = (v<0);

        %calculation of interaction terms for turning
        rmr_r = zeros(size(vp));
        lmr_r = zeros(size(vp));
        rml_r = zeros(size(vp));
        lml_r = zeros(size(vp));
        rmr_al = qal*vp*(gaussian(diffxp,sal,mal).*chi(2*sal-diffxp));
        lmr_al = qal*vm*(gaussian(diffxp,sal,mal).*chi(2*sal-diffxp));
        rml_al = qal*vp*(gaussian(diffxm,sal,mal).*chi(2*sal-diffxm));
        lml_al = qal*vm*(gaussian(diffxm,sal,mal).*chi(2*sal-diffxm));
        rmr_a = zeros(size(vp));
        lmr_a = zeros(size(vp));
        rml_a = zeros(size(vp));
        lml_a = zeros(size(vp));

        %create output matrix:
        output = vertcat(rmr_r,lmr_r,rml_r,lml_r,rmr_al,lmr_al, ...
                         rml_al,lml_al,-rmr_a,-lmr_a,-rml_a,-lml_a);

        signal_p = right_moving_flags*output;
        signal_m = left_moving_flags*output;

        %turning rates
        lambda_p = lambda1 + lambda2*f(normalization.*signal_p);
        lambda_m = lambda1 + lambda2*f(normalization.*signal_m);

        %update velocities
        rp = rand(1,N);
        rm = rand(1,N);
        v = chi(v).*sign(rp-lambda_p*dt)-chi(-v).*sign(rm-lambda_m* ...
                                                       dt);

        %%%%%%%%%%%%%%% Moving half-STEP

        %get new directions
        vp = (v>0);
        vm = (v<0);

        %calculation of interaction terms for moving
        rmr_r = qr*vp*(gaussian(diffxp,sr,mr).*chi(2*sr-diffxp));
        lmr_r = qr*vm*(gaussian(diffxp,sr,mr).*chi(2*sr-diffxp));
        rml_r = qr*vp*(gaussian(diffxm,sr,mr).*chi(2*sr-diffxm));
        lml_r = qr*vm*(gaussian(diffxm,sr,mr).*chi(2*sr-diffxm));
        rmr_al = zeros(size(vp));
        lmr_al = zeros(size(vp));
        rml_al = zeros(size(vp));
        lml_al = zeros(size(vp));
        rmr_a = qa.*vp*(gaussian(diffxp,sa,ma).*chi(2*sa-diffxp));
        lmr_a = qa.*vm*(gaussian(diffxp,sa,ma).*chi(2*sa-diffxp));
        rml_a = qa.*vp*(gaussian(diffxm,sa,ma).*chi(2*sa-diffxm));
        lml_a = qa.*vm*(gaussian(diffxm,sa,ma).*chi(2*sa-diffxm));

        %speed output (attraction - repulsion + 0*alignment)
        speedoutput = vertcat(-rmr_r,-lmr_r,-rml_r,-lml_r,rmr_al, ...
                              lmr_al,rml_al,lml_al,rmr_a, lmr_a,rml_a,lml_a);

        speed_signal_p = right_moving_flags*speedoutput;
        speed_signal_m = left_moving_flags*speedoutput;

        %speeds
        gamma_p = gamma*g(normalization.*speed_signal_p);
        gamma_m = gamma*g(normalization.*speed_signal_m);

        %update positions
        x = x + gamma_p.*vp.*dt - gamma_m.*vm.*dt;
        %periodic boundaries
        x = mod(x,L);
        %store positions to plot later
        traj = [traj; x];

    end
    y = traj;
end
