function [profil_rgb,xcenter,ycenter,radius] = image_analysis_2_rot(im,angle,k,coef)
% find center and radius of the interference pattern (presumed to be
% circular). Also returns the horizontal profile on which will be computed
% the thickness (profil_rgb).
% - im: image to be analyzed
% - angle: the profile is measured on the horizontal line passing through 
% the center of the pattern. If you want to measure the profile on a 
% specific line, you need to rotate the image with an angle specified in 
% degrees. 
% - k: number of lines on which the intensity is averaged
% - coef: coefficient used to binarize the image and find the center of the
% interference patter. Need to be adjusted for each case. In my case:
% coef=1.1 

% coef = 1.2;

im2 = im;%rgb2hsv(im);
BW = im2bw(im2(:,:,1),coef*graythresh(im2(:,:,1)));
BW = imfill(BW,'holes');
stats = regionprops('table',BW,'Centroid','Area', 'MajorAxisLength','MinorAxisLength');
stats2 = stats(stats.Area == max(stats.Area),:);
xcenter = stats2.Centroid(:,1);
ycenter = stats2.Centroid(:,2);

diameter = mean([stats2.MajorAxisLength stats2.MinorAxisLength],2);
radius = diameter/2;

% figure; 
% imshow(im)
% hold on ;
% plot(xcenter,ycenter,'o')
% viscircles(stats2.Centroid,radius);

% subplot(3,1,2)
% hold on ;
% % on moyenne sur plusieurs lignes et on vérifie que cela ne change pas
% % grand chose
% for k = 1:2%5:30;
im = rotateAround(im,ycenter,xcenter,angle); % vérifier ordre x,y

figure; 
imshow(im)
hold on ;
% plot(xcenter,ycenter,'o')
% viscircles(stats2.Centroid,radius);

% k = 3; 
% profil2 = mean(im(round(ycenter)-k:round(ycenter+k),round(xcenter-radius):round(xcenter+radius),:));
% % profil2 = im(round(ycenter),:,:);
% profil2 = double(profil2)/255;
% profil_rgb = squeeze(profil2);
if k == 0
    profil2 = im(round(ycenter):round(ycenter),round(xcenter-radius):round(xcenter+radius),:);
    profil2 = double(profil2)/255;
    profil_rgb = squeeze(profil2);
    profil_rgb = imadjust(profil_rgb,[min(profil_rgb(:)) min(profil_rgb(:)) min(profil_rgb(:)) ; max(profil_rgb(:)) max(profil_rgb(:)) max(profil_rgb(:))]);
else 
    profil2 = mean(im(round(ycenter)-k:round(ycenter)+k,round(xcenter-radius):round(xcenter+radius),:));
    profil2 = double(profil2)/255;
    profil_rgb = squeeze(profil2);
    profil_rgb = imadjust(profil_rgb,[min(profil_rgb(:)) min(profil_rgb(:)) min(profil_rgb(:)) ; max(profil_rgb(:)) max(profil_rgb(:)) max(profil_rgb(:))]);
% plot(profil)
% end
% 
% subplot(3,1,3)
% plot(profil,'.');
% minp = min(profil_rgb(:)) ;
% maxp = max(profil_rgb(:)) ;

end