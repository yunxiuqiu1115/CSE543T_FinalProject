function [xs, ys, raceinfos] = buildTrainCellArrays(data, years, states)
    xs = cell(1000,1);
    ys = cell(1000,1);
    raceinfos = cell(1000,1);
    counter = 1;
    for i = 1:numel(years)
       for j = 1:numel(states)
          [x, y, candidateNames, v, pvi, experienced] = getRaceCandidateData(data, years(i), states(j));
          if isempty(x), continue; end
          for k = 1:numel(x)
             xs{counter} = x{k};
             ys{counter} = y{k};
             raceinfos{counter} = {years(i), states(j), candidateNames(k), v(k), pvi(k), experienced(k)};
             counter = counter + 1;
          end
       end
    end

    counter = counter - 1;
    xs = xs(1:counter);
    ys = ys(1:counter);
    raceinfos = raceinfos(1:counter);
end