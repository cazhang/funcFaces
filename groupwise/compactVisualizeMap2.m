% visualize maps of a specific pair, but with different method
function [h1, h2, h3, h4, h5, h6] = compactVisualizeMap2(mesh1, basis1, mesh2, basis2, Fmap)
% Visualizes the map from the given vertex.
nMethod = length(Fmap.maps);
dest = cell(nMethod, 1);
fmap = cell(nMethod, 1);

for m = 1:length(Fmap.maps)    
    if isfield(Fmap,'ppMaps')
        fmap{m} = Fmap.ppMaps{m};
        dest{m} = Fmap.ppDestinations{m};
    else
        fmap{m} = Fmap.maps{m};
        dest{m} = Fmap.destinations{m};
    end
    
        %set(h,'WindowStyle','docked','name',getTitle(Fmap,m));    
end

[h, h1,h2, h3, h4, h5, h6] = compactShowFaceMap(mesh1, basis1, mesh2, basis2, fmap, dest);
set(h,'WindowStyle','normal', 'name', getTitle(Fmap));  


function titleStr = getTitle(Fmap)
%titleStr = Fmap.method{1};
titleStr = 'Eval:';
for m=1:length(Fmap.maps)
    if isfield(Fmap,'params') && ~isempty(Fmap.params)
        titleStr = sprintf('%s [(%g)', titleStr, Fmap.params{m});
    end
    if isfield(Fmap,'ppMaps')
        titleStr = sprintf('%s + %s', titleStr, Fmap.ppMethod);
        if isfield(Fmap,'ppParams') 
            titleStr = sprintf('%s (%g)', titleStr, Fmap.ppParams{1});
        end        
        error = Fmap.ppErrors(m);
    else
        error = Fmap.errors(m);
    end
    titleStr = sprintf('%s, Error: %g]',titleStr, error);
end

