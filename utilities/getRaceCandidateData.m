function [xs, ys, candidateNames, vs, pvis, experienceds]= getRaceCandidateData(data, year, state)
    oneRaceData = data(data.cycle==year & strcmp(data.state, state), :);
    if ~size(oneRaceData,1), xs = []; ys = []; candidateNames = []; vs = []; pvis = []; experienceds = []; return ; end
    oneRaceData(:, ["cycle", "state", "pollster", "Democrat"]) = [];
    
    % divide data into candidate group
    [G, candidates_identifiers] = findgroups(oneRaceData.Candidateidentifier);
    oneRaceData(:, "Candidateidentifier") = [];
    for i=1:numel(candidates_identifiers)
        candidate_identifier = candidates_identifiers{i};
        % remove year/state
        candidate_identifier = extractAfter(candidate_identifier,6);
        candidate_identifier(~ismember(candidate_identifier, ['A':'Z', 'a':'z', '0':'9'])) = '';
        allCandidatesData.(candidate_identifier) = sortrows(oneRaceData(G==i, :), ["daysLeft","samplesize"]);
    end
    
    candidateNames = fieldnames(allCandidatesData);
    nc = numel(candidateNames);
    xs = cell(nc,1);
    ys = cell(nc,1);
    vs = zeros(nc, 1);
    pvis = zeros(nc, 1);
    experienceds = zeros(nc, 1);
    for i=1:nc
        % define training data
        ts = allCandidatesData.(candidateNames{i}).daysLeft;             % daysLeft
        ns = allCandidatesData.(candidateNames{i}).numberSupport;        % numberSupport
        ss = allCandidatesData.(candidateNames{i}).samplesize;           % samplesize 
        ds = allCandidatesData.(candidateNames{i}).Republican;           % Republican or not
        ps = ns./ss;                                                     % polling proportions
        is = allCandidatesData.(candidateNames{i}).pollsteridx;          % pollster indexes
        vs(i) = allCandidatesData.(candidateNames{i}).Percentage_of_Vote_won_x(1); % votes won
        xs{i} = [ts, ps, ss, is, ds];
        % candidateid = candidateid + 1;
        ys{i} = ps;
        pvis(i) = allCandidatesData.(candidateNames{i}).pvi(1);
        experienceds(i) = allCandidatesData.(candidateNames{i}).experienced(1);
    end
end