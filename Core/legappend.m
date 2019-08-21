function [legend_h,object_h,plot_h,text_strings] = legappend(newStrings,varargin)
%LEGAPPEND appends new entries to the end of a legend by deleting the
%current legend and recreating a new, similar legend. 
% 
%% Syntax
% 
% legappend('new legend entry') 
% legappend('new entry 1','new entry 2',...,'new entry N') 
% legappend('') 
% legappend('','',...,'')
% [legend_h,object_h,plot_h,text_strings] = legappend(...)
% 
%% Description 
% 
% legappend('new legend entry') appends an existing legend with "new
% legend entry".
% 
% legappend('new entry 1','new entry 2',...,'new entry N') adds several
% new entries to the legend. 
% 
% legappend('') deletes the last entry from the legend. 
% 
% legappend('','',...,'') deletes the last several entries from the
% legend.
% 
% [legend_h,object_h,plot_h,text_strings] = legappend(...) returns legend_h, the
% handle of the new legend; object_h, handles of the line, patch, and
% text graphics objects used in the legend; plot_h, handles of the lines
% and other objects used in the plot; and text_strings, a cell array of
% the text strings used in the legend. Note that for new legend entries,
% legappend does not add entries to a current legend, but deletes the
% current legend and recreates a new one. As a result, the legend handle
% will change with each new-entry use of legappend.  The legend handle
% does not change when legappend is used to delete an entry. 
% 
% 
%% Author Info
% This function was created by Chad A. Greene of the Institute for
% Geophysics, The University of Texas at Austin, July 2014. 
% 
% See also legend.

h =  findobj(gcf,'Type','axes','Tag','legend');

prop.boxon = get(h,'visible');
prop.loc = get(h,'location'); 
prop.color = get(h,'color'); 
prop.orient = get(h,'Orientation'); 



allDatah = flipud(get(gca,'children')); 
str = get(h,'String'); 

if exist('varargin','var') 
    newStrings = [newStrings,varargin];
end
deleteEntries = sum(cellfun('isempty',newStrings));
if isempty(newStrings) 
    deleteEntries = 1; 
end

if ~deleteEntries
    if iscell(newStrings)
        for k = 1:length(newStrings) 
            str{end+1}=newStrings{k}; 
        end
    end 
    if ~iscell(newStrings)
        str{end+1}=newStrings; 
    end


    [legend_h,object_h,plot_h,text_strings] = legend(h,allDatah,str);

    if strcmpi({prop.boxon},'off')
        legend boxoff
    end

    set(legend_h,'location',prop.loc,'color',prop.color,'Orientation',prop.orient)


end

if deleteEntries
    set(h,'String',str(1:end-nargin))
    [legend_h,object_h,plot_h,text_strings] = legend;
end


if nargout==0
    clear legend_h object_h plot_h text_strings
end


end

