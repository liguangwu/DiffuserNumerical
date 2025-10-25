function [coolingtime, Temp_final, isbreak] = rate2time(app,coolingrate,coolingpath,u0,T0,DiffCoef,dx,f_dt,tol_percent)
%The diffusion continues when the temperature decreases until T reaches a
%threshold value capable of freezing the diffusion profile
isbreak=0;
% tol_percent=1e-8; %tolerance in percent
R=8.314;
D0=exp(DiffCoef{1}); E=DiffCoef{3};
coolingtime = NaN(length(coolingrate),1); %final time when the diffuison is frozen
Temp_final = NaN(length(coolingrate),1);
% wb=uiprogressdlg(app.DiffuserGUI,'Title','Calculating corresponding timescales...','Indeterminate','on','Cancelable', 'on');
wb=uiprogressdlg(app.DiffuserGUI,'Title','Calculating corresponding timescales...','Indeterminate','on');

parfor cr=1:length(coolingrate)
    % if wb.CancelRequested
    %     isbreak=1;
    %     break
    % end
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

        if T<=273.15
            T=Temp(end-2);
            t=time(end-2);
            Temp(end-1:end)=[];
            time(end-1:end)=[];
            u_fit(:,end-1:end)=[];
            u=u_fit(:,end);
            Di(:,end)=[];
            D=Di(:,end);
            dt=dt0; %decreases T more slowly
            dt_fix=1;
            continue
        end

        D=ones(numel(u),1)*D0*exp(-E*1000/R/T);
        Di = cat(2,Di,D);

        if ~dt_fix
            dt=(dx)^2/(2*max(D))*f_dt; %Cooling makes D becomes smaller and dt larger, which may result in next T<0
        else
            dt=dt*1.2;
        end

        %difference between last and this step
        if i>1
            delta=sum((u-u_fit(:,end-1)).^2);
            tol=sum((u_fit(:,end-1)*tol_percent).^2); %tolerance
            if delta<=tol %|| T<(800+273.15)
                isstop = 1; %when to stop
                coolingtime(cr) = t;
                Temp_final(cr) = T;
            end
        end

    end
end

close(wb);