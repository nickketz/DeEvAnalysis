function [dep,Eout] = deev_dependency(res,pair,c)

%%function to calculate behavioural dependency measure
% Aidan J Horner 04/2015

% input:
% res       =   MxN matrix of M 'events' and N 'retrieval trials'
% pair      =   pair of retrieval trials to calculate dependency 
%               (e.g., cue location retrieve object and cue location retrieve person)
%               example = [1 2] - uses first two columns to build contigency table for analysis
% optional inputs:
% c         =   number of choices (i.e., c-alternative forced choice) - for
%               estimating level of guessing
%               default = 6
%
% output:
% dep       =   dependency measure for [data independent_model dependent_model dependent_model+guess]
% Eout      =   episodic factor for each event

%% housekeeping

if nargin   < 3
    c       = 6;                                            % set 'c' to 6 if not defined by user
end

%% calculate dependency for data
 
res2        = res(:,pair);                                  % create column Mx2 matrix for retrieval trials defined by 'pair'
dep(1)      = sum(sum(res2,2)~=1)/size(res2,1);         	% calculate dependency for data

%% calculate dependency for independent model

acc         = mean(res2,1);                                 % calculate accuracy for each retrieval type
dep(2)      = ((acc(1)*acc(2))+((1-acc(1))*(1-acc(2))));    % calculate dependenct for independent model

%% calculate dependency for dependent model
for guess = [0 1]
    Eout        = nan(size(res2,1),1);
    cont        = nan(size(res2,1),2,2);                        % create matrix for dependent model probabilities
    g           = (1-nanmean(res(:)))*(c/(c-1));                % calculate level of guessing !!!modified mean to be nanmean for my res structure includes nans!!!
    b           = mean(res); b(:,pair) = nan;                   % calculate average performance
    for i       = 1:size(res2,1)                                % loop through all event
        a       = res(i,:); a(:,pair) = nan;                    % calculate event specific performance
        Eout(i) = nanmean(a)/nanmean(b);                        % calculate ratio of event / average performance (episodic factor) --nk: shouldn't this be ratio of sums? I guess its equivalent if they have the same number of elements
        E = Eout(i);
        for p = 1:2;
            if E*acc(p)>1
                P(p) = 1;
            else
                if guess == 1
                    P(p) = (E*(acc(p)-(g/c)))+(g/c);
                elseif guess == 0
                    P(p) = E*acc(p);
                end
            end
        end
        cont(i,1,1) = P(1)*P(2);
        cont(i,1,2) = (1-P(1))*P(2);
        cont(i,2,1) = P(1)*(1-P(2));
        cont(i,2,2) = (1-P(1))*(1-P(2));
    end
    cont2       = squeeze(sum(cont));                           % create contingency table
    %dep(3)      = (cont2(1,1)+cont2(2,2))/sum(cont2(:));        % calculate dependency for dependent model
    tmp         = (cont2(1,1)+cont2(2,2))/sum(cont2(:));        % calculate dependency for dependent model
    if guess
        dep(4) = tmp;
    else
        dep(3) = tmp;
    end
end
