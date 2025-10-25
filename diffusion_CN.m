function [time, misfit, rd, x_grid, u_fit, u_fit_conv, isbreak, Ibd]= diffusion_CN(...
    profile_x,profile_C,Weight,initial_x,initial_C,Nx,f_dt,T,DiffCoef,deconv_parameters)
% isothermal conditions, calculate timescale
% Crank-Nicolson solution for for D=D0*exp(-E/RT)
isbreak=0;
%set parameters------------------------------------------------------------
D0=exp(DiffCoef{1}); E=DiffCoef{3}; %kJ/mol
R = 8.314; %gas constant
Lx=max(profile_x)-min(profile_x); %diffusion profile length
dx=Lx/Nx; %grid size in x
x_grid=min(profile_x)+(0:Nx)*dx; %x values on the grid
x_grid=x_grid';
%initial composition
u = interp1(initial_x,initial_C,x_grid,"linear","extrap");
D=ones(numel(u),1)*D0*exp(-E*1000/R/T);
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
    
    a0=[-4*D(2:end-1)+D(3:end)-D(1:end-2);0;0];
    b0=[1;8*dx^2/dt+8*D(2:end-1);1];
    c0=[0;0;-4*D(2:end-1)-D(3:end)+D(1:end-2);];
    a1=-a0;
    b1=[1;8*dx^2/dt-8*D(2:end-1);1];
    c1=-c0;
    A=spdiags([a0 b0 c0],-1:1,length(u),length(u));
    B=spdiags([a1 b1 c1],-1:1,length(u),length(u));
    u=A\(B*u);
    u_fit=cat(2,u_fit,u);

    % D=ones(numel(u),1)*D0*exp(-E*1000/R/T);
    % dt=(dx)^2/(2*max(D))*f_dt;
    % drawnow
    % hf1=figure(1);
    % pos=hf1.Position;
    % plot(profile_x,profile_C,'o');
    % hold("on")
    % plot(x_grid,u,'r-');
    % hold("off")
      
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
    profile_pred = interp1(x_grid,u2,profile_x,"nearest","extrap");
    misfit = cat(1, misfit, sum(Weight.*(profile_pred-profile_C).^2) );
    rd = cat(2,rd,profile_pred-profile_C);
    % hf2=figure(2);
    % hf2.Position=[pos(1)+pos(3) pos(2) pos(3) pos(4)];
    % plot(time, misfit)

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

% [~,time_step]=min(misfit)