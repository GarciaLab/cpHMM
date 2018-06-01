function ChangeColorPBoC2(PlotHandle,ColorProperty)

CurrentColor=get(PlotHandle,ColorProperty);

if isempty(strmatch(CurrentColor,'none'))&isempty(strmatch(CurrentColor,'auto'))
    if sum(CurrentColor==[1,0,0])==3
        set(PlotHandle,ColorProperty,[213,108,85]/255)
    elseif sum(CurrentColor==[0,0,1])==3
        set(PlotHandle,ColorProperty,[115,142,193]/255)
    elseif sum(CurrentColor==[0,1,0])==3 
        set(PlotHandle,ColorProperty,[122,169,116]/255)
    elseif sum(CurrentColor==[1,1,0])==3 %Yellow
        set(PlotHandle,ColorProperty,[234,194,100]/255)
    elseif sum(CurrentColor==[0,1,1])==3 %Cyan
        set(PlotHandle,ColorProperty,[108,188,233]/255)
    elseif sum(CurrentColor==[1,0,1])==3 %Magenta
        set(PlotHandle,ColorProperty,[208,109,171]/255) 
    end
end