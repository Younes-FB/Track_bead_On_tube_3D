%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain1=2; %%% channel 1 intensity gain 
gain2=4; %%% channel 2 intensity gain
smfactor=3; % smothing factor 
pausetime=0.01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentfolder=pwd;
[file,Folder] = uigetfile('*.tif');
fullFileName = [Folder file];
outputname=file(1:end-4);
      prompt = {'Slices(z)','Frame(t)','Frame interval (s)'};
      dlg_title = 'Parameters:';
      num_lines = 1;
      defaultans = {'28','25','13'};
      answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
      Slices=str2num(answer{1});
      Frame_nr=str2num(answer{2});
      Fstep=str2num(answer{2});
xb=ones;yb=ones;zb=ones;xf=ones;yf=ones;zf=ones;          
hA=figure; set(hA,'Position',[100 100 1000 500]); 
for j = 1:Frame_nr
    clear im0 Imtube Imbead imyzb imyzt I1 I2 I3 I4 
    for i = 1:Slices*2
         n=Slices*2*(j-1)+i;
         im0(:,:,i) = imread(fullFileName,'Index',n);
    end
     Imtube=im0(:,:,1:2:end-1);
     Imbead=im0(:,:,2:2:end);
     imyzt=permute(Imtube,[1 3 2]);
     imyzb=permute(Imbead,[1 3 2]); 
     [sy,sx,sz]=size(Imtube);
%%%%%%%%%%%%%%%%%%%%% find the center of the bead
         windowWidth = 5; % Whatever.  Bigger number for smoother images.
         kernel = ones(windowWidth) / windowWidth^2.5;
         I1 = conv2(double(max(Imbead,[],3)), kernel, 'same');
         I2 = conv2(double(max(Imtube,[],3)), kernel, 'same');
         [yb0,xb0]=find(I1==max(max(I1)));  
         xb(j)=round(mean(xb0));yb(j)=round(mean(yb0));
      subplot(2,2,1);
      C1 = imfuse(I1.*gain2,I2.*gain1,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
      imshow(C1);
        set(gca,'xaxisLocation','top')
        xlabel('x (pixel)');ylabel('y (pixel)');
        hold on
        plot(xb(j),yb(j),'o','MarkerFaceColor','w')      
%%%%%%%%%%%%%%%%%%%%%%%% find location of center of the bead and filopodium
    xstart=round(xb0)-2;xend=round(xb0)+2;
    if xstart<1
        xstart=1;
    elseif xend>sx
        xend=sx;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         I3 = conv2(double(max(imyzb(:,:,xstart:xend),[],3)), kernel, 'same');
         I4 = conv2(double(max(imyzt(:,:,xstart:xend),[],3)), kernel, 'same');
         [yb00,zb0]=find(I3==max(max(I3)));
         [yt0,zt0]=find(I4==max(max(I4)));
  yf(j)=round(mean(yt0));zf(j)=round(mean(zt0));
  zb(j)=round(mean(zb0));yb(j)=round(mean(yb00));
  subplot(2,2,2);
  C2 = imfuse(I3.*gain2,I4.*gain1,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
         imshow(C2);
         set(gca,'xaxisLocation','top')
        xlabel('z (pixel)');ylabel('y (pixel)');
     hold on
        plot(zb(j),yb(j),'O','MarkerFaceColor','w','Markersize',8) 
     hold on
       plot(zf(j),yf(j),'O','MarkerFaceColor','yellow','Markersize',4)  
       pause(pausetime)
end
       %%%% flip z direction(if image is not already flipped in imageJ, active this part)
         %        zb=-zb+max(zb)+min(zb);zf=-zf+max(zb)+min(zb);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct the tube movement
       dx=xb;dy=yf-yf(1);dz=zf-zf(1);
       xbn=xb;ybn=yb-dy;zbn=zb-dz;
       yfn=yf-dy;zfn=zf-dz;xfn=xbn;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zbn = smooth(zbn,smfactor);ybn = smooth(ybn,smfactor);xbn = smooth(xbn,smfactor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 3D trajectory
subplot(2,2,3);
plot3(xfn,yfn,zfn,'.-g','MarkerSize',30,'linewidth',10);
          xlabel('x (pixel)');xlim([min(xbn)-5 max(xbn)+5]);
          ylabel('y (pixel)');ylim([min(ybn)-5 max(ybn)+5]);
          zlabel('z (pixel)');zlim([min(zbn)-5 max(zbn)+5]);
 c = 1:numel(xbn);      
 hold on
 h = surface([xbn(:), xbn(:)], [ybn(:), ybn(:)],[zbn(:), zbn(:)],...
    [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none','linewidth',1);
colormap('jet')
hold on 
T=0:1:length(xbn)-1;
scatter3(xbn,ybn,zbn, 30, T'*Fstep, 'filled')
view(-15,35)
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 2D view from the tip side of the filamen
subplot(2,2,4);
plot(zfn,yfn,'.g','MarkerSize',30);
hold on
scatter(zbn, ybn, 20, T'*Fstep, 'filled')    
          xlabel('z (pixel)');xlim([min(zbn)-5 max(zbn)+5]);
          ylabel('y (pixel)');ylim([min(ybn)-5 max(ybn)+5]);
colormap('jet')
cl=colorbar;
cl.Label.String = 'Time(s)';
box on
grid on
u = gradient(zbn);
v = gradient(ybn);
scale = 0.6;
quiver(zbn, ybn,u,v,scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% save data
    cd(Folder)
    Namef=[file(1:end-4) '_track_Plots'];       
    print('-dpng',Namef)
  
column_names = {'  x_tube', '  y_tube ', 'z_tube  ','x_bead  ','y_bead ','z_bead  '};
Res=[xfn;yfn;zfn;xbn';ybn';zbn'];
fname = sprintf([file(1:end-4) '_positions.txt']);
fileID = fopen(fname,'w');
fprintf(fileID, '%6s ', column_names{:});
fprintf(fileID,'\n %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f',Res);
fclose(fileID);

            