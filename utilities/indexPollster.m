function [data,pollster2idx] = indexPollster(data, threshold)
% Index pollster above threshold
%
% Note: pollster below threshold will be indexed into one category
%       this category is indexed by max(pollsteridx)

    % selecting pollsters with sufficient polls
    if isa(threshold,'double')
        pollsters = cell2table(tabulate(data.pollster), 'VariableNames', {'pollster', 'count','percent'});
        pollsters.percent = [];
        selectedpollsters = pollsters(pollsters.count >= threshold, :);

        pollster2idx = containers.Map;
        for i=1:size(selectedpollsters,1)
            pollster2idx(selectedpollsters.pollster{i}) = i;
        end

        nfirm = size(selectedpollsters,1) + 1;
        for i=1:size(pollsters,1)
            if ~isKey(pollster2idx, pollsters.pollster{i}), pollster2idx(pollsters.pollster{i}) = nfirm; end
        end

        idxs = zeros(size(data,1),1);
        for i=1:size(data,1)
           idxs(i) = pollster2idx(data.pollster{i}); 
        end

        data.pollsteridx = idxs;
    %     writetable(data, file_path);
    else
        pollster2idx = threshold;

        idxs = zeros(size(data,1),1);
        for i=1:size(data,1)
            if ~isKey(pollster2idx, data.pollster{i})
                idxs(i) = 41;
            else
                idxs(i) = pollster2idx(data.pollster{i});
            end
        end

        data.pollsteridx = idxs;
    end
end