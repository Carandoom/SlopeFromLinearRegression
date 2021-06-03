%	
%	Need dot separated numbers in the original txt file
%   
%	To add: converter from comma to dot numbers
%   Take the DELTA from the beginning and end of the slope
%   Adjust the slope to not take into account the flat parts
%   Add filter for time section to work in
%   
%   DATA STRUCTURE:
%       Text files is imported, the first column is stored into x and used
%       as x values whereas all the other columns are stored into T and
%       used as y values.
%   
%   Example on a negative slop to find, upper graph is the trace to fit
%   and the bottom on is the first derivative. The red dots are the ones
%   taken into account for the linear regression.
%       Youtube link: https://youtu.be/a8gxQElWiVE
%   
%   Christopher HENRY - V1 - 2021-02-11
%   

%%  To add: 
%   Add an error handling for mdl generation if no data
%   Add the delta measurement
%   Find a way to take data from previous macro ?
%   add a data selection tool ?
%   find a way to rearange the data and generate conditions to group the
%   dishes in
%
%   /!\ Add a way to convert comma data into dot data /!\


function FindSlope()

close all
clear
clc

[T,x,path,filename] = OpenFile();
[f,ax,ylimPlot] = CreateFigure(T,x);

waitfor(findobj('Tag','bNxt')); %   Pause until range selected
tic

[dydx] = FindPositiveSlope(f,T,x);

%  Extract longest slope and fit linear regression
Parameters = NaN(2,size(T,2));
ax(2).Visible = 'on';
for i = 1:size(T,2)
    if isnan(T(1,i))
        continue
    end
    try
        pause(1);
        plot(ax(1),x,T(:,i),'Color', 'b');
        plot(ax(2),x,dydx(:,i),'Color', 'b');
        ylimPlot = ylim(ax(1));
        [Parameters] = FindLongestSlope(ax,T,x,dydx,i,Parameters);
        ylim(ax(1),[ylimPlot(1) ylimPlot(2)]);
    catch
    end
        disp([int2str(i) ' out of ' int2str(size(T,2))])
end
ylim([ylimPlot(1) ylimPlot(2)]);

%  Plot Parameters
boxplot(ax(2),Parameters(1,:));
hold on
scatter(ones(1,size(T,2)),Parameters(1,:));
hold off

% Write data in txt
fileID = fopen([path 'Slope_' filename],'w');
fprintf(fileID,'%f\n',Parameters(1,:));
fclose(fileID);

toc
end


function [T,x,path,filename] = OpenFile()
    cd 'C:\Users\henryc\Desktop\GitHub\Matlab find slope'
    [filename, path] = uigetfile('*.*', ...
        'MultiSelect', 'off');
    T = readtable([path filename]);
    T = table2array(T);
    x = T(:,1);
    T = T(:,2:end);
    [T] = MedSmoothData(T);
%     [T] = SmoothData(T);
end
function [T] = MedSmoothData(T)
    MedLen = 5;
    T = medfilt2(T,[MedLen 1]);
end
function [T] = SmoothData(T)
    for i = 1:size(T,2)
        T(:,i) = smooth(T(:,i));
    end
end

function [f,ax,ylimPlot] = CreateFigure(T,x)
    f = figure();
    ax(1) = subplot(2,1,1);
    ax(2) = subplot(2,1,2);
    ax(2).Visible = 'off';
    
    [ylimPlot] = PlotAllData(ax,T,x);
    AddSlider(f,x(1),x(end));
    CreateCheckPosNegSlope(f);
    CreateValidationButton(f);
end
function [ylimPlot] = PlotAllData(ax,T,x)
    plot(ax(1),...
        x,...
        T,...
        'Color', 'b');
    ylimPlot = ylim(ax(1));
    xlim(ax(1),[x(1), x(end)]);
end
function [] = AddSlider(f,ValMin,ValMax)
    txt1 = uicontrol(f, ...
        'Style', 'text', ...
        'Tag', 'txt1', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'String', ValMin, ...
        'Units', 'normalized', ...
        'Position', [0.02, 0.3966, 0.1, 0.06] ...
        );
    uit1 = uicontrol(f, ...
        'Style', 'slider', ...
        'Tag', 'uit1', ...
        'Min', ValMin, ...
        'Max', ValMax, ...
        'Value', ValMin, ...
        'Units', 'normalized', ...
        'Position', [0.13, 0.4, 0.775, 0.06] ...
        );
    addlistener(uit1,...
        'ContinuousValueChange',...
        @(hObject,event) ChangeVal(hObject,event,uit1.Value,txt1));
    txt2 = uicontrol(f, ...
        'Style', 'text', ...
        'Tag', 'txt2', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'String', ValMax, ...
        'Units', 'normalized', ...
        'Position', [0.02, 0.2966, 0.1, 0.06] ...
        );
    uit2 = uicontrol(f, ...
        'Style', 'slider', ...
        'Tag', 'uit2', ...
        'Min', ValMin, ...
        'Max', ValMax, ...
        'Value', ValMax, ...
        'Units', 'normalized', ...
        'Position', [0.13, 0.3, 0.775, 0.06] ...
        );
    addlistener(uit2,...
        'ContinuousValueChange',...
        @(hObject,event) ChangeVal(hObject,event,uit2.Value,txt2));
end
function [] = CreateCheckPosNegSlope(f)
    CheckBox1 = uicontrol(f, ...
        'Style', 'checkbox', ...
        'Tag', 'CheckBox1', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'String', 'Negative Slope', ...
        'Units', 'normalized', ...
        'Position', [0.42, 0.185, 0.3, 0.08] ...
        );
end
function [] = CreateValidationButton(f)
    bNxt = uicontrol(f,...
        'Style', 'pushbutton', ...
        'String', 'Validate Range', ...
        'FontSize',12, ...
        'Tag', 'bNxt', ...
        'Units', 'normalized', ...
        'Position', [0.36, 0.1, 0.3, 0.08], ...
        'CallBack', @(hObject,event) DeleteElements(f) ...
        );
end
function DeleteElements(f)
    uit1 = findobj('Tag','uit1');
    uit2 = findobj('Tag','uit2');
    txt1 = findobj('Tag','txt1');
    txt2 = findobj('Tag','txt2');
    bNxt = findobj('Tag','bNxt');
    CheckBox1 = findobj('Tag','CheckBox1');
    if round(uit1.Value)>round(uit2.Value)
        warndlg('Min value must be lower than Max value','Error range value')
        return
    end
    f.UserData = [round(uit1.Value) round(uit2.Value) CheckBox1.Value;];
    uit1.delete
    uit2.delete
    txt1.delete
    txt2.delete
    bNxt.delete
    CheckBox1.delete
end

function [] = ChangeVal(hObject,event,NewVal,ObjText)
    ObjText.String = int2str(NewVal);
end

function [dydx] = FindPositiveSlope(f,T,x)
    [FirstN, LastN] = GetRangeSlope(f,x);
    [dydx] = FirstDerivativeSmooth(T,x);
    if f.UserData(3) == 1
        dydx = -dydx;
    end
    dydx(dydx<=0) = 0;
    dydx([1:FirstN LastN:end],:) = 0;
end
function [FirstN, LastN] = GetRangeSlope(f,x)
    FirstN = f.UserData(1);
    LastN = f.UserData(2);
    [~, FirstN] = min(abs(x-FirstN));
    [~, LastN] = min(abs(x-LastN));
end
function [dydx] = FirstDerivativeSmooth(T,x)
    dydx = zeros(size(T,1), size(T,2));
    for i = 1:size(T,2)
        dydx(:,i) = gradient(T(:,i)) ./ gradient(x);
        dydx(:,i) = smooth(dydx(:,i),10);
    end
end

function [Parameters] = FindLongestSlope(ax,T,x,dydx,i,Parameters)
    zpos = find(~[0 dydx(:,i)' 0]);
    [~, grpidx] = max(diff(zpos));
    
    [Parameters] = LinearFit(ax,T,x,dydx,zpos,grpidx,i,Parameters);
end
function [Parameters] = LinearFit(ax,T,x,dydx,zpos,grpidx,i,Parameters)
    val1 = zpos(grpidx);
    val2 = zpos(grpidx+1)-2;
    [FitValY,ThrGaussFit] = HalfNormGaussFDeriv(x,dydx,i,val1,val2);
    
    mdl = fitlm(x(find(FitValY>=ThrGaussFit)+val1),...
        T(find(FitValY>=ThrGaussFit)+val1,i));
    
    Parameters(1,i) = mdl.Coefficients.Estimate(2);
    PlotLinearFit(ax,T,x,dydx,val1,val2,i,mdl,FitValY,ThrGaussFit)
end
function [FitValY,ThrGaussFit] = HalfNormGaussFDeriv(x,dydx,i,val1,val2)
    NormFit = fit(x(val1:val2),dydx(val1:val2,i),'gauss2');
    FitValY = NormFit(x(val1:val2));
    ThrGaussFit = min(FitValY)+(max(FitValY)-min(FitValY))/2;   
end
function [] = PlotLinearFit(ax,T,x,dydx,val1,val2,i,mdl,FitValY,ThrGaussFit)
    line(ax(1),x(find(FitValY>=ThrGaussFit)+val1),...
        T(find(FitValY>=ThrGaussFit)+val1,i),...
        'Marker','o',...
        'Color','r',...
        'LineStyle','none');
    line(ax(1),x(val1:val2),...
        mdl.Coefficients.Estimate(2).*x(val1:val2)+mdl.Coefficients.Estimate(1),...
        'Color','k');
    
    line(ax(2),x(find(FitValY>=ThrGaussFit)+val1),...
        dydx(find(FitValY>=ThrGaussFit)+val1,i),...
        'Marker','o',...
        'Color','r',...
        'LineStyle','none');
end


%%  DEBUG

function DEBUG()
    cd 'C:\Users\henryc\Desktop\GitHub\Matlab find slope'
    T = readtable('DataTest2.txt');
    T = table2array(T);
    x = T(:,1);
    T = T(:,2:end);
    
    test = fit(x(val1:val2),dydx(val1:val2,i),'gauss2');
    plot(test(x(val1:val2)));

end
