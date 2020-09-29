function [xs, ys, candidateNames, vs, pvis, experienceds, parties]= getRaceCandidateData(data, year, state)

% Obtain polling/fundemental features and election race metadata of given year/state
% Results in the form of cell array (length equals to number of candidate in this year/state race)
% input: 
%   - data: dataset that contains experience,
%           cycle,state,Candidateidentifier,Percentage_of_Vote_won_x,
%           samplesize,daysLeft,numberSupport,pollster,Republican,Democrat
%   - year: election year of to-be-collected features and metadata 
%   - state: election state of to-be-collected features and metadata
% 
% output:
%   - xs: cell array of polling/fundemental features.
%          Each element contains days before election, polling proportion and samplesize
%   - ys: cell array of pollings. Each element contains polling proportion.
%   - raceinfos: cell array of election race metadata. Each element
%   - candidateNames: names of candidates
%   - vs: actual vote shares
%   - pvis: cook partisan voting indices
%   - experience: whether candidates has served in office
%   - parties: indicators of candidate parties (-1 if republican, 1 if democratic, 0 if third party)

    % obtain data given year/state
    oneRaceData = data(data.cycle==year & strcmp(data.state, state), :);
    
    % return empty arrays in there is no election for that year/state combination
    if ~size(oneRaceData,1), xs = []; ys = []; candidateNames = []; vs = []; pvis = []; experienceds = []; parties = []; return ; end
    oneRaceData(:, ["cycle", "state", "pollster"]) = [];
    
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
    
    % initialize output arrays
    candidateNames = fieldnames(allCandidatesData);
    nc = numel(candidateNames);
    xs = cell(nc,1);
    ys = cell(nc,1);
    vs = zeros(nc, 1);
    parties = zeros(nc, 1);
    pvis = zeros(nc, 1);
    experienceds = zeros(nc, 1);
    for i=1:nc
        ts = allCandidatesData.(candidateNames{i}).daysLeft;             % daysLeft
        ns = allCandidatesData.(candidateNames{i}).numberSupport;        % numberSupport
        ss = allCandidatesData.(candidateNames{i}).samplesize;           % samplesize 
        rs = allCandidatesData.(candidateNames{i}).Republican;           % Republican or not
        ds = allCandidatesData.(candidateNames{i}).Democrat;             % Democrat or not
        ps = ns./ss;                                                     % polling proportions
        vs(i) = allCandidatesData.(candidateNames{i}).Percentage_of_Vote_won_x(1); % actual vote shares
        % ds, rs is the same across one candidate
        % indicators of candidate parties (-1 if republican, 1 if democratic, 0 if third party)
        if ds(1)==1
            parties(i) = 1;
        elseif rs(1) == 1
            parties(i) = -1;
        else
            parties(i) = 0;
        end
            
        xs{i} = [ts, ps, ss];
        ys{i} = ps;
        pvis(i) = allCandidatesData.(candidateNames{i}).pvi(1);
        experienceds(i) = allCandidatesData.(candidateNames{i}).experienced(1);
    end
end