function [coolingrate, coolingtime, Temp_final, misfit_final, rd_final, x_grid, u_final, u_final_conv, isbreak, Ibd] ...
    = diffusion_CN_Pl_rate(profile_x,profile_C,Weight,initial_x,initial_C,XAn,Nx,f_dt,T,coolingpath,DiffCoef,tol_percent,fO2,P,aSiO2,deconv_parameters)
% non-isothermal conditions, calculate cooling rates and timescales
% Crank-Nicolson solution for D=f(T,P,fO2,C...)
isbreak=0;
% tol_percent=1e-8; %tolerance in percent
%try isothermal to find best fit and boundary*************************************************
%***********************************************************************************************
[time, misfit, ~, ~, ~, ~, ~, Ibd]= diffusion_CN_Pl(...
    profile_x,profile_C,Weight,initial_x,initial_C,XAn,Nx,f_dt,T,DiffCoef,fO2,P,aSiO2,deconv_parameters);
misfit_95cl=max(misfit(Ibd)); %if misfit is lower than this value, we consider it as best fit
%cooling rate cannot be larger than this value
switch coolingpath
    case 'linear cooling'
        dTdtmax=(T-273.15)/time(Ibd(1)); %T=T-dT_dt*dt;
    case 'exponential cooling'
        dTdtmax=log(T/273.15)/time(Ibd(1)); %T=T*exp(-dT_dt*dt);
    case 'parabolic cooling'
        dTdtmax=(T-273.15)/time(Ibd(1))^2; %T=T-dT_dt*dt^2;
end
%***********************************************************************************************

%set parameters----------------------------------------------------------------
Lx=max(initial_x)-min(initial_x); %diffusion profile length
dx=Lx/Nx; %grid size in x
x_grid=min(initial_x)+(0:Nx)*dx; %x values on the grid
x_grid=x_grid';
%initial composition
u = interp1(initial_x,initial_C,x_grid,"linear","extrap");
XAn_grid = interp1(initial_x,XAn,x_grid,"linear","extrap");
u0 = u; T0 = T;
%find suitable values of dTdt---------------------------------------------------------
dTdt = dTdtmax; %max dTdt calculated by max T/min timescale
isstop=0; %find the 95c.l. of  dTdt
counts=0; %modeling of dTdt will make the diffusion curve move across two 95c.l. boundaries of best fit
coolingrate=[]; %store suitable dTdt
coolingtime=[]; %store corresponding timescale
Temp_final=[]; %store final timescale
u_final = [];
u_final_conv=[]; %for convolution
misfit_final=[];
rd_final = [];

i=0;
step=0.8; %dTdt(i+1)=dTdt(i)*step
backtrial=1;
Ibd=[NaN NaN]; %95c.l. of best fit
while ~isstop
    i=i+1;
    [iswithbd, timescale, Ci, Ci_conv, Ti, misfiti, rd] = goodrate_Pl(dTdt, coolingpath, ...
        u0, XAn_grid, T0, DiffCoef, fO2, P, aSiO2, dx, x_grid, profile_x, profile_C, Weight, misfit_95cl, f_dt, tol_percent, deconv_parameters);
    
    coolingrate=cat(1, coolingrate, dTdt);
    coolingtime=cat(1, coolingtime, timescale);
    u_final=cat(2, u_final, Ci);
    u_final_conv=cat(2, u_final_conv, Ci_conv);
    Temp_final=cat(1,Temp_final,Ti);
    misfit_final=cat(1, misfit_final, misfiti);
    rd_final=cat(1, rd_final, rd);

    if iswithbd %first reach the 95c.l. of best fit
        if counts==0
            counts=1;
            Ibd(1)=i; %misfit move from high rate to low rate
        end
        if backtrial
            %go back two steps
            if length(coolingrate)>1
                backstep=2;
                coolingrate(end-backstep+1:end)=[];
                coolingtime(end-backstep+1:end)=[];
                u_final(:,end-backstep+1:end)=[];
                u_final_conv(:,end-backstep+1:end)=[];
                misfit_final(end-backstep+1:end)=[];
            else 
                backstep=1;
                coolingrate=[];
                coolingtime=[];
                u_final=[];
                u_final_conv=[];
                misfit_final=[];
            end
            dTdt=dTdt/step.^(backstep)/0.95; %be sure that the new step is the same as below
            i=i-1-backstep;
            %make step smaller within 95c.l. of best fit
            step=0.95;
            backtrial=0;
            counts=0;
        end
    end
    if ~iswithbd && counts==1 %first beyond the 95c.l. of best fit
        Ibd(2)=i-1; %misfit move from high rate to low rate
        % counts=2;
        % isstop=1;
        break
    end
    dTdt=dTdt*step;
end


% [~,timestep]=min(misfit_final)
