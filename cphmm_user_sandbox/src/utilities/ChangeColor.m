function ChangeColor(PlotHandle,ColorProperty)

CurrentColor=get(PlotHandle,ColorProperty);

if isempty(strmatch(CurrentColor,'none'))
    if sum(CurrentColor==[1,0,0])==3
        set(PlotHandle,ColorProperty,[236,32,36]/255)
    elseif sum(CurrentColor==[0,0,1])==3
        set(PlotHandle,ColorProperty,[57,83,163]/255)
    elseif sum(CurrentColor==[0,1,0])==3 
        set(PlotHandle,ColorProperty,[106,189,70]/255)
    elseif sum(CurrentColor==[1,1,0])==3 %Yellow
        set(PlotHandle,ColorProperty,[255,242,0]/255)
    elseif sum(CurrentColor==[0,1,1])==3 %Cyan
        set(PlotHandle,ColorProperty,[50,153,204]/255)
    elseif sum(CurrentColor==[1,0,1])==3 %Magenta
        set(PlotHandle,ColorProperty,[237,38,144]/255) 
    end
end