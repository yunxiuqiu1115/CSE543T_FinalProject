function [xs, ys] = getCellArrayData(data, years, states)
    allx = cell(numel(years),numel(states));
    ally = cell(numel(years),numel(states));
    counter = 1;
    for i=1:numel(years)
       for j=1:numel(states)
            oneRaceData = data(data.cycle==years(i) & strcmp(data.state, states{j}), :);
            if ~size(oneRaceData,1), continue; end
            oneRaceData(:, ["cycle", "state", "pollster"]) = [];

            % divide data into candidate group
            [G, candidates_identifiers] = findgroups(oneRaceData.Candidateidentifier);
            oneRaceData(:, "Candidateidentifier") = [];
            for k=1:numel(candidates_identifiers)
                candidate_identifier = candidates_identifiers{k};
                % remove year/state
                candidate_identifier = extractAfter(candidate_identifier,6);
                candidate_identifier(~ismember(candidate_identifier, ['A':'Z', 'a':'z', '0':'9'])) = '';
                allCandidatesData.(candidate_identifier) = sortrows(oneRaceData(G==k, :), ["daysLeft","samplesize"]);
            end

            % building models for each candidate
            candidateNames = fieldnames(allCandidatesData);
            allx{i,j} = cell(numel(candidateNames), 1);
            ally{i,j} = cell(numel(candidateNames), 1);
            for k=1:numel(candidateNames)
                % define training data
                ts = allCandidatesData.(candidateNames{k}).daysLeft;             % daysLeft
                ns = allCandidatesData.(candidateNames{k}).numberSupport;        % numberSupport
                ss = allCandidatesData.(candidateNames{k}).samplesize;           % samplesize 
                ps = ns./ss;                                                     % polling proportions
                is = allCandidatesData.(candidateNames{k}).pollsteridx;          % pollster indexes
                x = [ts, ps, ns, is];
                y = ps;
                allx{i,j,k} = x;
                ally{i,j,k} = y;
                counter = counter + 1;
            end
       end
    end
   
    counter = counter - 1;
    xs = cell(counter, 1);
    ys = cell(counter, 1);
    
    count = 1;
    for i=1:numel(allx)
       if isempty(allx{i}), continue; end
       xs{count} = allx{i};
       ys{count} = ally{i};
       count = count + 1;
    end
    count = count - 1;
    assert(count == counter);
end