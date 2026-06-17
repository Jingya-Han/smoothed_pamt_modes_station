function [ppdf_kde, pamt_kde] = KDE(x, data, bandwidth, dbin)
% kernel density estimation in a faster speed
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Remove NaNs and zeroes from data
    logdata = log(data(:));
    valid_indices = ~isnan(logdata);   % Find valid data indices
    logdata = logdata(valid_indices);  % Remove NaNs from logdata
    data = data(valid_indices);        % Remove NaNs and corresponding values from data

    n = length(data);  % Number of data points
    m = length(x);     % Number of query points

    % Precompute log(x)
    x_log = log(x);

    % Create n x m matrices of log(x) and logdata
    x_log_matrix = repmat(x_log, n, 1);       % Repeat x_log as rows (n x m)
    logdata_matrix = repmat(logdata, 1, m);    % Repeat logdata as columns (n x m)

    % Vectorize kernel calculation using matrix operations
    kernel_matrix = exp(-0.5 * ((x_log_matrix - logdata_matrix) / bandwidth).^2) / (sqrt(2 * pi) * bandwidth);  % (n x m)

    % Frequency distribution
    f_hat = sum(kernel_matrix, 1) / n;  % Sum over columns, normalized by n (1 x m)
    ppdf_kde = f_hat * dbin;
    ppdf_kde(1) = 1 - sum(ppdf_kde(2:end));  % Ensure normalization of the first bin

    % Amount distribution
    kernel_matrix_p = kernel_matrix .* data;  % Element-wise multiplication (n x m)
    p_hat = sum(kernel_matrix_p, 1) / n;      % Sum over columns, normalized by n (1 x m)
    p_hat(1) =0; % set wet-day precip as 0
    
    % Renormalize the amount distribution by daily mean precipitation
    meanp = nanmean(data(:)); 
    pamt_kde = p_hat * (meanp / sum(p_hat));  % Normalize

    
end
