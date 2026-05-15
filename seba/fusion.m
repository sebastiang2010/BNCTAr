clc 

Dosis3=Dosis_boro; 

Dosis4=permute(Dosis3,[2 3 1]);






% for i=1:88 
%        fusionada_1(:,:,i)=imfuse(A(:,:,i),B(:,:,i),'montage');
% end 
% 
% figure(200)
% imshow(fusionada_1(:,:,57),[]);
% 
% % 
% % [optimizer, metric] = imregconfig('monomodal');
% % % Alinear las imágenes
% % %  tform = imregtform(A, B,'rigid',optimizer,metric);
% % % % 
% % % % % Aplicar transformación afín a la imagen de dosis
% % %  img_dosis_alineada = imwarp(A, tform, 'OutputView', imref3d(size(B)));
% % 
% % %if isfield(tform,'T')
% %     %dosis_aligned = imwarp(A, tform);
% %     % Asignar pesos a las imágenes antes de fusionarlas
% %     peso_dosis = 0.8;
% %     peso_ct = 0.2;
% %     
% % 
% %     fusionada=zeros(109,330,88);
% %     % Fusionar las imágenes con ponderación de intensidad
% %     for i=1:88 
% %        fusionada(:,:,i)=imfuse(img_dosis_alineada(:,:,i),B(:,:,i),'montage');
% %     end 
% %     % Mostrar la imagen resultante
% %     figure(100)
% %     imshow(fusionada(:,:,57),[]);
% % %else
% % %    disp('No se pudo obtener una transformación válida');
% % %end
% 
% 
% A1=A(end:-1:1,:,:);
% 
% for i=1:88 
%     fusionada_3(:,:,i)=imfuse(A1(:,:,i),B(:,:,i),'montage');
% end 
%  figure(300)
%  for i=1:88
%    imshow(fusionada_3(:,:,i),[]);
%  end 
% 
% 
% dose=A1(:,:,57);
% I=B(:,:,57);
% 
% I1=cat(3,I,I,I); 
% 
% figure(500)
% ax=gca;
% cla(ax)
% image(I1);
% hold on
% if patient.show>=1
%     % Generar un corte coronal
%     %h=imshow(dose_cor,[]);
%     %set(h,'AlphaData',0.4)
%     surf(dose,'FaceAlpha',0.4);
%     %set(patient.dose_surf,'ButtonDownFcn',@imgClick)
%     shading interp;
%     colormap jet
% end
% % txt=['Sagittal :  ',num2str(py)];
% % title(txt)
% axis square

