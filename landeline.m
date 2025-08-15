function [y] = landeline(x, sigF, sigM, optM)
    y = ((sigF^2)./(sigM^2) + 1)*(x - optM*((sigF^2)./(sigM^2)))
end