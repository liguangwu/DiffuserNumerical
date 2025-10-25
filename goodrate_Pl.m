function [iswithbd, t_final, C_final, C_final_conv, T_final, misfit_final, rd_final] = goodrate_Pl(coolingrate, coolingpath, ...
    u0, XAn_grid, T0, DiffCoef, fO2, P, aSiO2, dx, x_grid, profile_x, profile_C, Weight, misfit_95cl, f_dt, tol_percent, deconv_parameters)
%The diffusion continues when the temperature decreases until T reaches a
%threshold value capable of freezing the diffusion profile, then
%iswithbd: the final diffusion curve lies between the 95c.l. of best fit
%istooslow: cooling rate is too slow, the final diffusion curve lies beyond the best fit
%istoofast: cooling rate is too fast, the final diffusion curve hasn't reach the best fit
iswithbd = zeros(length(coolingrate), 1);
istooslow = zeros(length(coolingrate), 1);
istoofast = zeros(length(coolingrate), 1);
t_final = NaN(length(coolingrate),1); %final timescale when the diffuison is frozen
misfit_final = NaN(length(coolingrate),1); %final misfit when the diffuison is frozen
T_final = NaN(length(coolingrate),1); %final temperature when the diffusion ceases
C_final = NaN(length(u0),1); %final concentration profile when the diffusion ceases
C_final_conv = NaN(length(u0),1); %for convolution, final concentration profile when the diffusion ceases
rd_final = NaN(length(profile_C),1); %final residuals when the diffusion ceases

for cr=1:length(coolingrate)
    dT_dt = coolingrate(cr);
    %initial parameters for each cooling rate modeling
    u = u0; T = T0;
    isstop=0; %find the error curve
    i=0; %count time steps
    D=exp(DiffCoef(1)+DiffCoef(2)*log(fO2)+DiffCoef(3)*XAn_grid...
        +DiffCoef(4)/T+DiffCoef(5)*P+DiffCoef(6)*P/T+DiffCoef(7)*log(aSiO2));
    dt=(dx)^2/(2*max(D))*f_dt; %grid size of time, ensure stability
    dt_fix=0;
    dt0 = dt;
    time = [];
    Temp = [];
    misfit = [];
    Di=[];
    t = 0;
    u_fit=[];
    u_fit_conv=[]; %for convolution
    rd=[];
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
        %convolution=================
        deconvolution=deconv_parameters{strcmp(deconv_parameters(:,1), 'deconvolute'),2};
        if deconvolution
            u2 = conv_profile(x_grid,u,deconv_parameters);
        else
            u2 = u;
        end
        u_fit_conv = cat(2,u_fit_conv,u2);
        %=========================

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
                u_fit_conv(:,end-1:end)=[];
                Di(:,end)=[];
                misfit(end)=[];
                rd(:,end)=[];
                D=Di(:,end);
            elseif length(Temp)>1
                T=Temp(end-1);
                t=time(end-1);
                Temp(end)=[];
                time(end)=[];
                u_fit(:,end)=[];
                u_fit_conv(:,end)=[];
            end
            u=u_fit(:,end);
            u2=u_fit_conv(:,end);
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

        D=exp(DiffCoef(1)+DiffCoef(2)*log(fO2)+DiffCoef(3)*XAn_grid...
            +DiffCoef(4)/T+DiffCoef(5)*P+DiffCoef(6)*P/T+DiffCoef(7)*log(aSiO2));
        Di = cat(2,Di,D);

        if ~dt_fix
            dt=(dx)^2/(2*max(D))*f_dt; %Cooling makes D becomes smaller and dt larger, which may result in next T<0
        else
            dt=dt*1.2;
        end

        %misfit calculation
        profile_pred = interp1(x_grid,u2,profile_x,'nearest');
        misfit = cat(1, misfit, sum(Weight.*(profile_pred-profile_C).^2) );
        rd = cat(2, rd, profile_pred-profile_C);

        %difference between last and this step
        if i>1
            delta=sum((u2-u_fit_conv(:,end-1)).^2);
            tol=sum((u_fit_conv(:,end-1)*tol_percent).^2); %tolerance
            if delta<=tol %|| T<(800+273.15)
                isstop = 1; %when to stop
                t_final(cr) = t;
                misfit_final(cr) = misfit(end);
                T_final(cr) = T;
                C_final(:,cr) = u;
                C_final_conv(:,cr) = u2;
                rd_final(:,cr) = rd(:,end);
                if misfit(end)<=misfit_95cl %when the diffusion ceases, will the profile lcoate between the 95c.l. of best fit
                    iswithbd(cr)=1;
                    if ~isunix %not available on web app
                        figure(1);
                        hold on
                        x=Temp-273.15;
                        plot(x,misfit,'-','Color',[0.8 0.8 0.8]); %show closure temeprature
                        plot(x(end),misfit(end),'k.', 'MarkerSize', 10)
                        hold off
                        xlabel('Temperature (\circC)')
                        ylabel('RSS changes with T for each suitable cooling rate')
                        title(['Check whether msifit is constant when T drops to the end,'...
                            newline, 'i.e., whether diffusion profile is frozen.', newline,...
                            'If not, users have to lower the tolerance that defines when misfit is constant'])
                        % x=time;
                        % plot(x,misfit,'-','Color',[0.8 0.8 0.8]); %show closure temeprature
                        % [changePoints, ~] = findchangepts(misfit, 'Statistic', 'mean', 'MinThreshold', 0.1);
                        % if ~isempty(changePoints)
                        %     stable_start = changePoints(end); %stable since last change point
                        %     t_final(cr) = time(stable_start);
                        %     misfit_final(cr) = misfit(stable_start);
                        %     T_final(cr) = Temp(stable_start);
                        %     C_final(:,cr) = u_fit(:,stable_start);
                        %     C_final_conv(:,cr) = u_fit_conv(:,stable_start);
                        %     rd_final(:,cr) = rd(:,stable_start);
                        %     plot(x(stable_start), misfit(stable_start), 'k.', 'MarkerSize', 10);
                        % end
                        % hold off
                        % xlabel('Time (s)')
                        % ylabel(['RSS changes with time', newline, ' for each suitable cooling rate (each line)'])
                        % title(['Check whether msifit is constant when time elapses,'...
                        %     newline, 'i.e., whether diffusion profile is frozen.', ...
                        %     newline, 'Black spots show the frozen time'])
                    end
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