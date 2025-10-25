function [time, misfit, rd, x_grid, u_fit, u_fit_conv, isbreak, Ibd]= diffusion_CN_Ol(...
    profile_x,profile_C,Weight,initial_x,initial_C,Nx,f_dt,T,DiffCoef,fO2,P,aSiO2,deconv_parameters)
% isothermal conditions, calculate timescale
% Crank-Nicolson solution for D=f(T,P,fO2,C...)
isbreak=0;
%set parameters------------------------------------------------------------
Lx=max(initial_x)-min(initial_x); %diffusion profile length
dx=Lx/Nx; %grid size in x
x_grid=min(initial_x)+(0:Nx)*dx; %x values on the grid
x_grid=x_grid';
%initial composition
u = interp1(initial_x,initial_C,x_grid,"linear","extrap");
D=exp(DiffCoef(1)+DiffCoef(2)*log(fO2)+DiffCoef(3)*u...
    +DiffCoef(4)/T+DiffCoef(5)*P+DiffCoef(6)*P/T+DiffCoef(7)*log(aSiO2));
dt=(dx)^2/(2*max(D))*f_dt; %grid size of time, ensure stability
%Prepare misfit calculation between data and model:
t = 0;
time = [];
u_fit=[];
u_fit_conv=[]; %for convolution
misfit = [];
rd=[]; %residuals

%% finite difference modeling
isstop=0; %find the error curve
i=0; %count time steps
x=[x_grid; flipud(x_grid); x_grid(1)];
Ibd=[NaN NaN]; %95c.l. of best fit
while ~isstop
    i=i+1;
    t = t+dt; % update time (in sec)
    time = cat(1,time,t); % store time
    
    a0=-4*D(2:end-1)+D(3:end)-D(1:end-2);
    b0=8*dx^2/dt+8*D(2:end-1);
    c0=-4*D(2:end-1)-D(3:end)+D(1:end-2);
    d0=u(1:end-2).*(-a0)+u(2:end-1).*(8*dx^2/dt-8*D(2:end-1))+u(3:end).*(-c0);
    %a*C(i-1) + b*C(i) + c*C(i+1) = d
    a=[0; a0; 0];
    b=[1; b0; 1];
    c=[0; c0; 0];
    d=[u(1); d0; u(end)];
    %solve the tridiagonal matrix
    u=tridiag(a,b,c,d);
    u_fit=cat(2,u_fit,u);

    D=exp(DiffCoef(1)+DiffCoef(2)*log(fO2)+DiffCoef(3)*u...
        +DiffCoef(4)/T+DiffCoef(5)*P+DiffCoef(6)*P/T+DiffCoef(7)*log(aSiO2));
    dt=(dx)^2/(2*max(D))*f_dt;
    
    %convolution=================
    deconvolution=deconv_parameters{strcmp(deconv_parameters(:,1), 'deconvolute'),2};
    if deconvolution
       u2 = conv_profile(x_grid,u,deconv_parameters);
    else
       u2 = u;
    end
    u_fit_conv = cat(2,u_fit_conv,u2);
    %=========================
    %misfit calculation
    profile_pred = interp1(x_grid,u2,profile_x,'nearest','extrap');
    misfit = cat(1, misfit, sum(Weight.*(profile_pred-profile_C).^2) );
    rd = cat(2,rd,profile_pred-profile_C);
    
    %extend steps
    if i>3
        pk=findpeaks(-misfit);
        if ~isempty(pk)
            %find the error curve
            [~,Imin]=min(misfit);
            y=[u_fit_conv(:,1); flipud(u_fit_conv(:,end)); u_fit_conv(1,1)];
            [in, on]=inpolygon(profile_x, profile_C, x, y);
            useful=in | on;

            for j=(Imin-1):-1:1
                Ibd(1)=j; %left bound site
                tmp=misfit; tmp(1:Imin)=[];
                [~,I]=min(abs(tmp-misfit(j))); %right bound site
                Ibd(2)=I+Imin;
                %how many spots are within the bound
                y=[u_fit_conv(:,j); flipud(u_fit_conv(:,Ibd(2))); u_fit_conv(1,j)];
                [in, on]=inpolygon(profile_x(useful), profile_C(useful), x, y);
                if (sum(in)+sum(on)) >= floor(0.95*sum(useful))
                    break
                end
            end

            if sum(useful)>1 && misfit(end)>misfit(j)
                isstop=1;
            end
        end
    end
end

% [~,timestep]=min(misfit)