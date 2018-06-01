function StandardFigure(PlotHandle,AxisHandle,varargin)

%Make it such that it can identify things like '-o'.

%Figure out the journal option in varargin
if isempty(varargin)
    Journal=0;
elseif strcmp(varargin{1},'PLoS')
    Journal=1;
elseif strcmp(varargin{1},'Inset')
    Journal=2;
elseif strcmp(varargin{1},'PLoSTrace')
    Journal=3;
elseif strcmp(varargin{1},'PLoSInset')
    Journal=4;
elseif strcmp(varargin{1},'Methods')
    Journal=5;
elseif strcmp(varargin{1},'MethodsInsert')
    Journal=6;
end


%Added this so that the log tick marks wouldn't go away.
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
%Keep the size the same
    
for i=1:length(PlotHandle)
    %if ~(get(PlotHandle(i),'LineStyle')=='none')
	
    
    %Line style
    try
        LineStyle=get(PlotHandle(i),'LineStyle');
        if (~isempty(strmatch(LineStyle,'-')))|...
                (~isempty(strmatch(LineStyle,'--')))|...
                (~isempty(strmatch(LineStyle,'-.')))|...
                (~isempty(strmatch(LineStyle,'-o')))
            set(PlotHandle(i),'LineWidth',1)
        end
    end


    

    %Marker style
    try
        if get(PlotHandle(i),'Marker')=='.'
            set(PlotHandle(i),'MarkerSize',15)
            set(PlotHandle(i),'LineWidth',1)
        elseif (get(PlotHandle(i),'Marker')=='d')|...
            (strcmp(get(PlotHandle(i),'Marker'),'square'))|...
            (get(PlotHandle(i),'Marker')=='o')
            set(PlotHandle(i),'MarkerSize',5)
            set(PlotHandle(i),'LineWidth',1)
            if ~isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'none'))
                set(PlotHandle(i),'MarkerFaceColor','w')
            end
        end
    end

    %Color section
  
    if strcmp(get(PlotHandle(i),'Type'),'hggroup')
        try ChangeColor(PlotHandle(i),'FaceColor')
        end
        try ChangeColor(PlotHandle(i),'Color')
        end
    elseif strcmp(get(PlotHandle(i),'Type'),'errorbarseries')
        
    else
        %Why do I have this auto thingy?
        if isempty(strmatch(get(PlotHandle(i),'Color'),'auto'))
            ChangeColor(PlotHandle(i),'Color')
        end

        if isempty(strmatch(get(PlotHandle(i),'MarkerEdgeColor'),'auto'))
            ChangeColor(PlotHandle(i),'MarkerEdgeColor')
        end

        if isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'auto'))
            ChangeColor(PlotHandle(i),'MarkerFaceColor')
        end
    end

%     
%     
%     for j=1:length(ColorProperty)
%         
%         
%             
%         
%         [i,j]
%         
%         
%         CurrentColor=get(PlotHandle(i),ColorProperty{j});
%         
%         if strmatch(CurrentColor,'auto')
%             CurrentColor=get(PlotHandle(i),'MarkerFaceColor');
%         end
%            
%       
%     end
            
            
            
            
    
    
    
%     %Point style
%     if isempty(strmatch(get(PlotHandle(i),'LineStyle'),'none'))
%         set(PlotHandle(i),'LineWidth',1)
%     end
    
    %Marker style
%     if get(PlotHandle(i),'Marker')=='.'
%         set(PlotHandle(i),'MarkerSize',15)
%         set(PlotHandle(i),'LineWidth',1)
% 	elseif (get(PlotHandle(i),'Marker')=='d')|...
%             (get(PlotHandle(i),'Marker')=='s')|...
%             (get(PlotHandle(i),'Marker')=='o')
%         set(PlotHandle(i),'MarkerSize',5)
%         set(PlotHandle(i),'LineWidth',1)
%     end
end

if Journal==1
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',24,'FontName','Arial')
    set(get(AxisHandle,'YLabel'),'FontSize',24,'FontName','Arial')
    set(get(AxisHandle,'ZLabel'),'FontSize',24,'FontName','Arial')
    set(AxisHandle,'FontSize',20,'FontName','Arial')
elseif Journal==3
    set(gcf,'Position',[221   438   699   260*0.75])
    set(gca,'Position',[0.1300    0.1100    0.7750*0.85   0.8150*0.75])
    set(get(AxisHandle,'XLabel'),'FontSize',24,'FontName','Arial')
    set(get(AxisHandle,'YLabel'),'FontSize',24,'FontName','Arial')
    set(get(AxisHandle,'ZLabel'),'FontSize',24,'FontName','Arial')
    set(AxisHandle,'FontSize',20,'FontName','Arial')
elseif Journal==4
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',48,'FontName','Arial')
    set(get(AxisHandle,'YLabel'),'FontSize',48,'FontName','Arial')
    set(get(AxisHandle,'ZLabel'),'FontSize',48,'FontName','Arial')
    set(AxisHandle,'FontSize',40,'FontName','Arial')
elseif Journal==5;
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',20,'FontName','Arial')
    set(get(AxisHandle,'YLabel'),'FontSize',20,'FontName','Arial')
    set(get(AxisHandle,'ZLabel'),'FontSize',20,'FontName','Arial')
    set(AxisHandle,'FontSize',18,'FontName','Arial')
elseif Journal==6;
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',40,'FontName','Arial')
    set(get(AxisHandle,'YLabel'),'FontSize',40,'FontName','Arial')
    set(get(AxisHandle,'ZLabel'),'FontSize',40,'FontName','Arial')
    set(AxisHandle,'FontSize',36,'FontName','Arial')
elseif Journal==2
    set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',24,'FontName','Lucida Sans Regular')
    set(get(AxisHandle,'YLabel'),'FontSize',24,'FontName','Lucida Sans Regular')
    set(get(AxisHandle,'ZLabel'),'FontSize',24,'FontName','Lucida Sans Regular')
    set(AxisHandle,'FontSize',20,'FontName','Lucida Sans Regular')
 else
%     set(gca,'Position',[0.1300    0.1100    0.7750    0.8150])
    set(get(AxisHandle,'XLabel'),'FontSize',17,'FontName','Lucida Sans Regular')
    set(get(AxisHandle,'YLabel'),'FontSize',17,'FontName','Lucida Sans Regular')
    set(get(AxisHandle,'ZLabel'),'FontSize',17,'FontName','Lucida Sans Regular')
    set(AxisHandle,'FontSize',15,'FontName','Lucida Sans Regular')
% else
%     set(get(AxisHandle,'XLabel'),'FontSize',24,'FontName','Arial')
%     set(get(AxisHandle,'YLabel'),'FontSize',24,'FontName','Arial')
%     set(get(AxisHandle,'ZLabel'),'FontSize',24,'FontName','Arial')
%     set(AxisHandle,'FontSize',20,'FontName','Arial')
end
% set(AxisHandle,...
%     'Box','off',...
%     'Color',[0.8510,0.8510,0.8510],...
%     'XColor','w','YColor','w','ZColor','w',...
%     'TickLength',[0.02,0.05])


%Change the color of the axes
try ChangeColor(AxisHandle,'XColor')
end
try ChangeColor(AxisHandle,'YColor')
end
