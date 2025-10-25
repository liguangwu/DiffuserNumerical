function convoluted_downsampled = conv_profile(x,C,deconv_parameters)
%convolute the beam effect
%x=profile_x;
%C=profile_C;
%Example: profile_C_conv = conv(app,profile_x,profile_C);

convoluted_downsampled=[];
%prepare deconvolution-----------------------------------------
para_names = deconv_parameters(:,1);
beamshape=deconv_parameters{strcmp(para_names, 'beam shape'),2};
beamunit=deconv_parameters{strcmp(para_names, 'beam size unit'),2};
diameter=deconv_parameters{strcmp(para_names, 'diameter'),2};
width=deconv_parameters{strcmp(para_names, 'width'),2};
FWHM=deconv_parameters{strcmp(para_names, 'FWHM'),2};
FWHM_G=deconv_parameters{strcmp(para_names, 'FWHM Gauss'),2};
FWHM_L=deconv_parameters{strcmp(para_names, 'FWHM Lorentz'),2};

switch beamunit
    case 'nm'
        f=1e9;
    case 'um'
        f=1e6;
    case 'mm'
        f=1e3;
    case 'm'
        f=1;
end
try
    switch beamshape
        case 'circular'
            width=diameter/f; %circular, elliptical beam
            dx=min([width/4 ((max(x))-min(x))/(max([25 (max(size(x)))])*4)]);
            model_x=min(x)-(max(x)-min(x))/2:dx:max(x)+(max(x)-min(x))/2+dx;
            [X,Y] = meshgrid(0:dx:width,0:dx:width);
            centre=[mean(mean(X,1)) mean(mean(Y,2))];
            dist=sqrt(abs(X-centre(1)).^2+abs(Y-centre(2)).^2);
            dist(dist>width/2)=NaN;
            dist(~isnan(dist))=1;
            dist(isnan(dist))=0;
            proportions=sum(dist)./(sum(sum(dist)));
            proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
        case 'rectangle'
            width=width/f; %rectangle, square beam
            dx=min([width/4 ((max(x))-min(x))/(max([25 (max(size(x)))])*4)]);
            model_x=min(x)-(max(x)-min(x))/2:dx:max(x)+(max(x)-min(x))/2+dx;
            n_conv=round(width/mean(diff(model_x)));
            proportions=ones(1,n_conv);
            proportions=proportions./sum(proportions);
            proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
        case 'gauss'
            FWHM=FWHM/f; %full width at half maximum; Gaussian beam
            dx=min([FWHM/4 ((max(x))-min(x))/(max([25 (max(size(x)))])*4)]);
            model_x=min(x)-(max(x)-min(x))/2:dx:max(x)+(max(x)-min(x))/2+dx;
            sigma=FWHM/2.355;
            proportions=exp(-((model_x-mean(model_x)).^2./(2.*sigma.^2)));
            proportions=proportions./sum(proportions);
            proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
        case 'lorentz'
            FWHM=FWHM/f; %full width at half maximum; Lorentzian beam
            dx=min([FWHM/4 ((max(x))-min(x))/(max([25 (max(size(x)))])*4)]);
            model_x=min(x)-(max(x)-min(x))/2:dx:max(x)+(max(x)-min(x))/2+dx;
            proportions=(1./pi).*((0.5.*FWHM)./((model_x-mean(model_x)).^2+(0.5.*FWHM).^2));
            proportions=proportions./sum(proportions);
            proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
        case 'voigt'
            FWHM_L=FWHM_L/f; FWHM_G=FWHM_G/f; %full width at half maximum; Voigt beam
            FWHM=(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)); %get FWHM effective for voigt
            dx=min([FWHM/4 ((max(x))-min(x))/(max([25 (max(size(x)))])*4)]);
            model_x=min(x)-(max(x)-min(x))/2:dx:max(x)+(max(x)-min(x))/2+dx;
            proportions=((1.36603*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))-0.47719*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))^2+0.11116*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))^3)*((1./(pi*(0.5.*(FWHM_L+1e-10)))).*(((0.5.*(FWHM_L+1e-10))^2)./((model_x-mean(model_x)).^2+(0.5.*(FWHM_L+1e-10)).^2)))+(1-(1.36603*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))-0.47719*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))^2+0.11116*(FWHM_L/(0.5346*FWHM_L+sqrt(0.2166*FWHM_L^2+FWHM_G^2)))^3))*((1/(((FWHM_G+1e-10)/2.355)*sqrt(2*pi)))*exp(-0.5.*((model_x-mean(model_x))./((FWHM_G+1e-10)/2.355)).^2)));
            proportions=proportions./sum(proportions);
            proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
    end
catch
    errordlg('Please check the unit of beam size and distance','Warning');
end


original = interp1(x, C, model_x, 'linear', 'extrap');
% 直接使用conv，确保卷积核已归一化
proportions = proportions ./ sum(proportions);
convoluted = conv(original, proportions, 'same');
convoluted_downsampled = interp1(model_x,convoluted,x,'linear');
% plot(x,[C,convoluted_downsampled])
end