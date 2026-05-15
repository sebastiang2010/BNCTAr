function plotdose(img,dose,nimg,pmax)
global patient;

hold all

% dose(img==0)=NaN;

if(~isfield(patient,'hcont'))
    [C,patient.hcont] = contour(permute(100*(dose(:,:,nimg)/max(max(max(dose())))),[1 2 3]), 'LineColor', [0 0 0]);
    clabel(C,patient.hcont,'Color','white')
    patient.hcont.ContourZLevel = 200;
    set(patient.hcont,'ButtonDownFcn',@imgClick)
else
    
end

s1 = surf(permute(dose(:,:,nimg),[1 2 3]),'FaceAlpha',0.5);
set(s1,'ButtonDownFcn',@imgClick)
caxis([min(min(min(dose()))) max(max(max(dose())))]);
shading interp;

[m,n]=size(I);
levels=hist(reshape(I,[1 m*n]),256);
levels=find(levels)-1;

for i = 2:length(levels)
    BW=false(m,n);
    BW(permute(img(:,:,nimg),[1 2 3])==levels(i))=true;
    B=bwboundaries(BW,'noholes');
    for j = 1:length(B)
        bound = B{j};
        h = plot3(bound(:,2),bound(:,1),199*ones(length(bound),1),'k','LineWidth',1);
        set(h,'ButtonDownFcn',@imgClick)
    end
end

if (nimg==pmax(3))
    h = plot3(pmax(2),pmax(1),200,'k*');
    set(h,'ButtonDownFcn',@imgClick)
    text(pmax(2)+3,pmax(1)+3,200,'max','Color','k')
end
%hc = colorbar;
%set(get(hc,'title'),'string','Dose Rate [Gy/min]')
pbaspect([1 1 1])
colormap jet
axis off

hold off
end
