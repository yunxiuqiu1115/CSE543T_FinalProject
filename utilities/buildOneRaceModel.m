function [success, predPolls, trueVotes] = buildOneRaceModel(data, year, state, parms)
    [xs, ys, candidateNames, trueVotes] = getRaceCandidateData(data, year, state);
    if isempty(xs), success = -1; predPolls = []; trueVotes = []; return ; end
    
    % building models for each candidate
    nc = numel(candidateNames);
    predPolls = zeros(nc, 1);
    
    for i=1:nc
        % training model
        parms.i = i;
        [predPolls(i), fig] = oneCandidatePredict(xs{i}, ys{i}, parms);
        plot_title = year + " " + state + " " + candidateNames{i};
        title(plot_title);
        yearFolder = fullfile("plots", num2str(year));
        stateFolder = fullfile(yearFolder, state);
        if ~exist(yearFolder, 'dir')
            mkdir(yearFolder)
        end
        if ~exist(stateFolder, 'dir')
            mkdir(stateFolder)
        end
        filename = fullfile(stateFolder, plot_title + ".jpg");
        saveas(fig, filename);
        close;
        disp(plot_title + " predicted winning rate: " + predPolls(i));
        disp(plot_title + " actual votes won: " + trueVotes(i) + newline);
    end

    [~, p_idx] = max(predPolls);
    [~, t_idx] = max(trueVotes);

    if p_idx == t_idx
       disp("Successfully predicted!"); 
       success = 1;
    else
        success = 0;
        disp("Wrong predicton!"); 
    end

end
