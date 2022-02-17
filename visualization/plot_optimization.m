function optimal_params = plot_optimization(enerdata, paramdata, varargin)
%Plot Energy vs. params for a set of params
%enerdata is a vector of the energies at each SR bin
%paramdata is a matrix with each parameter in a column and same # of rows
%as enerdata
%varargin{1} is the name of the brewermap scheme to use.  default is
%'Spectral'


if nargin == 2
    colors = brewermap(size(paramdata,2),'Spectral');
elseif nargin == 3
    colors = brewermap(size(paramdata,2),varargin{1});
else
    error('Must use 2 or 3 arguments in plot_optimization') 
end

figure()
hold on
for v = 1:size(paramdata,2)
    if verLessThan('matlab', '9.11')
        for bin = 1:size(paramdata,1)
            scatter(paramdata(bin, v), enerdata(bin), 'MarkerFaceColor', colors(v,:), 'MarkerEdgeColor', colors(v,:), 'MarkerFaceAlpha', bin/size(enerdata,1));
        end
    else
        s = scatter(paramdata(:, v), enerdata, 'MarkerFaceColor', colors(v,:), 'MarkerEdgeColor', colors(v,:));
        s.AlphaData = (1:size(enerdata,1))'./size(enerdata,1);
        s.MarkerFaceAlpha = 'flat';
    end
end
xlabel('V param values')
ylabel('Energy')

optimal_params = paramdata(end,:);

end

