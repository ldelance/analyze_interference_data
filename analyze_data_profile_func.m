function [im_h] = analyze_data_profile_func(ImgDir,filename)
cd('C:\Users\Léa Delance\Documents\à sauvegarder\Manips\bike wheel\resultats juillet 2022\code_matlab')
addpath('C:\Users\Léa Delance\Documents\à sauvegarder\Manips\bike wheel\code_matlab_2303')
%%
% plage d'épaisseur pour le calcul théorique
h = [1:1:1000]; % ATTENTION, il ne faut pas changer le pas pour la définition de h

% calcul de l'indice optique
ncyOH = 1.453 ; % indice de réfraction du cyclopentanol
nC10 = 1.4097 ; % indice de réfraction du décane

n = 0.3 *nC10 + (1-0.3)*ncyOH ;

% calcul de la couleur théorique et conversion en HSV
plot_color = 0;
[th_color,th_hue,th_sat]=compute_color(h,n,plot_color);
th_color = squeeze(th_color);
%%
delta_idx = 80;
th_h = h;
% ImgDir = 'C:\Users\Léa Delance\Documents\à sauvegarder\Manips\bike wheel\résultats décembre 2022\last_image\300000cSt_800ppm\';
cd(ImgDir)
% listing = dir('*f*') ;

% for k =64%:length(listing)
im = imread(strcat(ImgDir,filename));
% im = imread(strcat(ImgDir,'221209_f1.tiff'));
% im_h = zeros(size(im,1),size(im,2));

[profil_rgb,xcenter,ycenter,radius] = image_analysis_2(im);

% radius = radius *1.06;

    profil2 = mean(im(round(ycenter)-3:round(ycenter)+3,round(xcenter-radius):round(xcenter+radius),:));
    profil2 = double(profil2)/255;
    profil_rgb = squeeze(profil2);
    profil_rgb = imadjust(profil_rgb,[min(profil_rgb(:)) min(profil_rgb(:)) min(profil_rgb(:)) ; max(profil_rgb(:)) max(profil_rgb(:)) max(profil_rgb(:))]);
    [prev_profil_idx_h] = calc_profil(profil_rgb,th_color,th_h, round(radius), delta_idx);


    im_h = th_h(prev_profil_idx_h);

% profil_rgb = imadjust(profil_rgb,[min(profil_rgb(:)) min(profil_rgb(:)) min(profil_rgb(:)) ; max(profil_rgb(:)) max(profil_rgb(:)) max(profil_rgb(:))]);
% [prev_profil_idx_h] = calc_profil(profil_rgb,th_color,th_h, round(radius), delta_idx);

% figure(101); hold on ;
% plot(im_h)

% end

function [prev_profil_idx_h] = calc_profil(profil_rgb,th_color,th_h, center, delta_idx)

prev_min = inf;
prev_profil_idx_h = zeros(length(profil_rgb),1);

for idx_h_init =1:length(th_h)
    
    [profil_idx_h] = calc_profil_with_hint(profil_rgb,center,idx_h_init,th_color,th_h,delta_idx);
    profil_h = th_color(profil_idx_h,:);
    di = color_distance_2(profil_h,profil_rgb);

    if di < prev_min
        prev_min = di;
        prev_profil_idx_h = profil_idx_h;
    end 
end

end



function [profil_idx_h] = calc_profil_with_hint(profil_rgb,idx_init,idx_h_init,th_color,th_h,delta_idx)

profil_idx_h = zeros(length(profil_rgb),1);
profil_idx_h(idx_init)=idx_h_init;

for k = idx_init+1:length(profil_rgb)


    idx_h_prev = profil_idx_h(k-1);
    idx_h_possible = NaN(length(th_h),1);
    
    if idx_h_prev - delta_idx >0 && idx_h_prev + delta_idx <= length(th_h)
    idx_h_possible(idx_h_prev - delta_idx : idx_h_prev + delta_idx) = 1 ; 
    elseif idx_h_prev - delta_idx <=0
    idx_h_possible(1 : idx_h_prev + delta_idx) = 1 ; 
    else
    idx_h_possible(idx_h_prev - delta_idx : length(th_h)) = 1 ;     
    end
    
    lsq = color_distance(profil_rgb(k,:),th_color);
    lsq = lsq .* idx_h_possible;
    [~,I] = min(lsq);
    profil_idx_h(k) = I;
    
end

for k = idx_init-1:-1:1


    idx_h_prev = profil_idx_h(k+1);
    idx_h_possible = NaN(length(th_h),1);
    
    if idx_h_prev - delta_idx >0 && idx_h_prev + delta_idx <= length(th_h)
    idx_h_possible(idx_h_prev - delta_idx : idx_h_prev + delta_idx) = 1 ; 
    elseif idx_h_prev - delta_idx <=0
    idx_h_possible(1 : idx_h_prev + delta_idx) = 1 ; 
    else
    idx_h_possible(idx_h_prev - delta_idx : length(th_h)) = 1 ;     
    end
    
    lsq = color_distance(profil_rgb(k,:),th_color);
    lsq = lsq .* idx_h_possible;
    [~,I] = min(lsq);
    profil_idx_h(k) = I;
    
end

end


function [lsq]=color_distance(profil_rgb,th_color)
    lsq = (profil_rgb(1)-th_color(:,1)).^2 + (profil_rgb(2)-th_color(:,2)).^2 +    (profil_rgb(3)-th_color(:,3)).^2;
% lsq = abs(profil_rgb(1)-th_color(:,1)) + abs(profil_rgb(2)-th_color(:,2)) +    abs(profil_rgb(3)-th_color(:,3));
end


function [di] = color_distance_2(profil_h,profil_rgb)
    d = profil_h-profil_rgb;
    di = sum(d(:,1).^2 + d(:,2).^2 + d(:,3).^2) ;

end

function [xunit,yunit] = coord_circle(x,y,r)
th = 0:pi/100:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end
end