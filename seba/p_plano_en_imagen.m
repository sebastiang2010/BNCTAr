% prueba cruces en imagenes 
clear 
clc 
close all 

load('mesh_2.mat')

I=patient.stack; 



x=10; %inicializo por si no lo toca 
y=20; 
z=80; 

I1=squeeze(I(x,:,:));
I2=squeeze(I(:,y,:));
I3=squeeze(I(:,:,z));


figure(1)
imshow(I1,[]);
hold on
set(gcf,'WindowButtonDownFcn',@f_imgClick)
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
plot([x x], [ylim(1) ylim(2)], 'r');
plot([xlim(1) xlim(2)], [y y], 'g');
scatter(x,y)

figure(2)
imshow(I2,[]);
hold on
set(gcf,'WindowButtonDownFcn',@f_imgClick)
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
plot([x x], [ylim(1) ylim(2)], 'r');
plot([xlim(1) xlim(2)], [y y], 'g');
scatter(x,y)

figure(3)
imshow(I3,[]);
hold on
set(gcf,'WindowButtonDownFcn',@f_imgClick)
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
plot([x x], [ylim(1) ylim(2)], 'r');
plot([xlim(1) xlim(2)], [y y], 'g');
scatter(x,y)


while true
    
    I1=squeeze(I(x,:,:));
    I2=squeeze(I(:,y,:));
    I3=squeeze(I(:,:,z));

    %analizar cual es z 

    figure(2)
    imshow(I1)
    hold on
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    plot([x x], [ylim(1) ylim(2)], 'r');
    plot([xlim(1) xlim(2)], [y y], 'g');
    scatter(x,y)
    set(gcf,'WindowButtonDownFcn',@f_imgClick)

    figure(3)
    imshow(I2)
    hold on
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    plot([x x], [ylim(1) ylim(2)], 'r');
    plot([xlim(1) xlim(2)], [y y], 'g');
    scatter(x,y)
    set(gcf,'WindowButtonDownFcn',@f_imgClick)

    figure(1)
    imshow(I3,[]);
    hold on
    set(gcf,'WindowButtonDownFcn',@f_imgClick)
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    plot([x x], [ylim(1) ylim(2)], 'r');
    plot([xlim(1) xlim(2)], [y y], 'g');
    scatter(x,y)


    % Espera a que se presione una tecla
    key = waitforbuttonpress;

    % Si la tecla presionada es "Esc", sal del bucle
    if strcmpi(get(gcf, 'CurrentKey'), 'escape')
        break;
    end

end
