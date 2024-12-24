function [iswithbd, t, C_final, T_final, misfit_final, rd_final] = goodrate(coolingrate, coolingpath, ...
    u0, T0, DiffCoef, dx, x_grid, profile_x, profile_C, Weight, misfit_95cl, f_dt, tol_percent)
%The diffusion continues when the temperature decreases until T reaches a
%threshold value capable of freezing the diffusion profile, then
%iswithbd: the final diffusion curve lies between the 95c.l. of best fit
%istooslow: cooling rate is too slow, the final diffusion curve lies beyond the best fit
%istoofast: cooling rate is too fast, the final diffusion curve hasn't reach the best fit
iswithbd = zeros(length(coolingrate), 1);
istooslow = zeros(length(coolingrate), 1);
istoofast = zeros(length(coolingrate), 1);
misfit_final = NaN(length(coolingrate),1); %final misfit when the diffuison is frozen
T_final = NaN(length(coolingrate),1); %final temperature when the diffusion ceases
C_final = NaN(length(u0),1); %final concentration profile when the diffusion ceases
rd_final = NaN(length(profile_C),1); %final residuals when the diffusion ceases
D0=exp(DiffCoef{1}); E=DiffCoef{3}; %kJ/mol
R = 8.314; %gas constant
for cr=1:length(coolingrate)
    dT_dt = coolingrate(cr);
    %initial parameters for each cooling rate modeling
    u = u0; T = T0;
    isstop=0; %find the error curve
    i=0; %count time steps
    D=ones(numel(u),1)*D0*exp(-E*1000/R/T);
    dt=(dx)^2/(2*max(D))*f_dt; %grid size of time, ensure stability
    dt_fix=0;
    dt0 = dt;
    time = [];
    Temp = [];
    misfit = [];
    Di=[];
    t = 0;
    u_fit=[];
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

        %update D
        switch coolingpath
            case 'linear cooling'
                T=T0-dT_dt*t;
            case 'exponential cooling'
                T=T0*exp(-dT_dt*t);
            case 'parabolic cooling'
                T=T0-dT_dt*t^2;
        end
        Temp = cat(1,Temp, T); %store temperature

        if T<=273.15 %go back to 1 or 2 steps
            if length(Temp)>2
                T=Temp(end-2);
                t=time(end-2);
                Temp(end-1:end)=[];
                time(end-1:end)=[];
                u_fit(:,end-1:end)=[];
                Di(:,end)=[];
                misfit(end)=[];
                rd(:,end)=[];
                D=Di(:,end);
            else
                T=Temp(end-1);
                t=time(end-1);
                Temp(end)=[];
                time(end)=[];
                u_fit(:,end)=[];
            end
            u=u_fit(:,end);
            dt=dt0; %decreases T more slowly
            dt_fix=1;
            continue
            % isstop = 1; %when to stop
            % misfit_final(cr) = misfit(end);
            % T_final(cr) = T;
            % C_final(:,cr) = u;
            % rd_final(:,cr) = rd;
            % return
        end

        D=ones(numel(u),1)*D0*exp(-E*1000/R/T);
        Di = cat(2,Di,D);

        if ~dt_fix
            dt=(dx)^2/(2*max(D))*f_dt; %Cooling makes D becomes smaller and dt larger, which may result in next T<0
        else
            dt=dt*1.2;
        end

        %misfit calculation
        profile_pred = interp1(x_grid,u,profile_x,'nearest');
        misfit = cat(1, misfit, sum(Weight.*(profile_pred-profile_C).^2) );
        rd = profile_pred-profile_C;

        %difference between last and this step
        if i>1
            delta=sum((u-u_fit(:,end-1)).^2);
            tol=sum((u_fit(:,end-1)*tol_percent).^2); %tolerance
            if delta<=tol
                isstop = 1; %when to stop
                misfit_final(cr) = misfit(end);
                T_final(cr) = T;
                C_final(:,cr) = u;
                rd_final(:,cr) = rd;
                if misfit(end)<=misfit_95cl %when the diffusion ceases, will the profile lcoate between the 95c.l. of best fit
                    iswithbd(cr)=1;
                    % figure(1);
                    % hold on
                    % plot(Temp-273.15,misfit) %show closure temeprature
                    % hold off
                    % xlabel('Temperature (\circC)')
                    % ylabel('misfit changes with T for each suitable cooling rate')
                    % title(['Check whether msifit is constant when T drops to the end,'...
                    %     newline, 'i.e., whether diffusion profile is frozen.', newline,...
                    %     'If not, users have to lower the tolerance that defines when misfit looks constant'])
                end
                if min(misfit)>=misfit(end) %whether reaches the best fit
                    istoofast(cr)=1;
                else
                    istooslow(cr)=1;
                end
            end
        end

    end
end