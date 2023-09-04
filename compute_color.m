function [th_color,th_hue,th_sat]=compute_color(h,n,plot_color)
% compute the theoretical observed depending on the film's thickness,
% refractive index and the spectrum of camera and lamp.
% h: thickness' range on which to make the computation
% n: refractive index of the liquid
% plot_color: =1 if you want to plot a figure with thickness to color
% conversion
%             =0 if you don't
% WARNING: This function requires to save previously the camera's sensor 
% sensitivity and lamp spectrum in separated .mat file. The files supplied
% with this function are adapted for a specific camera and lamp and are 
% only examples.



%%
% lg = 524;
% lb = 453;
% lr = 602;
% 
% dg = 0.5*(1-cos(4*pi*n*h/lg));
% db = 0.5*(1-cos(4*pi*n*h/lb));
% dr = 0.5*(1-cos(4*pi*n*h/lr));

%%

load('sensor sensitivity.mat');
load('lamp_spectrum.mat');


% 
    red = red(red(:,1)<740,:);
    green = green(green(:,1)<740,:);
    blue = blue(blue(:,1)<740,:);
    
%     red(:,2) = red(:,2)/max(red(:,2));
%     green(:,2) = green(:,2)/max(green(:,2));
%     blue(:,2) = blue(:,2)/max(blue(:,2));
%     
%     figure; plot(blue(:,1),blue(:,2)); hold on ;
%  plot(red(:,1),red(:,2))
% plot(green(:,1),green(:,2))
 lamp_red = interp1(lamp_spectrum(:,1),lamp_spectrum(:,2), red(:,1)); % need to check
lamp_green = interp1(lamp_spectrum(:,1),lamp_spectrum(:,2), green(:,1));
lamp_blue = interp1(lamp_spectrum(:,1),lamp_spectrum(:,2), blue(:,1));

  
    dr = zeros(1,length(h));
    dg = zeros(1,length(h));
    db = zeros(1,length(h));
        
    
for k = 1:length(h)
    
       dr(k) = trapz(red(:,1),lamp_red.*red(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./red(:,1))))/trapz(red(:,1),lamp_red.*red(:,2));
       dg(k) = trapz(green(:,1),lamp_green.*green(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./green(:,1))))/trapz(green(:,1),lamp_green.*green(:,2));
       db(k) = trapz(blue(:,1),lamp_blue.*blue(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./blue(:,1))))/trapz(blue(:,1),lamp_blue.*blue(:,2));
%        dr(k) = trapz(red(:,1),red(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./red(:,1))))/trapz(red(:,1),red(:,2)); % trapz is an integration function
%        dg(k) = trapz(green(:,1),green(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./green(:,1))))/trapz(green(:,1),green(:,2));
%        db(k) = trapz(blue(:,1),blue(:,2) .* 0.5.*(1-cos(4*pi*n*h(k)./blue(:,1))))/trapz(blue(:,1),blue(:,2));
end

%%
th_color(:,:,1) = dr;
th_color(:,:,2) = dg;
th_color(:,:,3) = db;

th_color_hsv = rgb2hsv(th_color);
th_hue = th_color_hsv(:,:,1);
th_sat = th_color_hsv(:,:,2);

if plot_color == 1
    im = zeros(1,length(h),3);
for k = 1:length(h)
   
    im(1,k,1) = dr(k);
    im(1,k,2) = dg(k);
    im(1,k,3) = db(k);

    
end

im2 = [im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im; im];

figure; imshow(im2);
axis on
xticks(h(100:100:end))
set(gca,'ytick',[])
set(gca,'yticklabel',[])
end
end