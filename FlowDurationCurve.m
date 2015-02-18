% Flow Duration Curve - Solomon Vimal
clear all; 
close all;
data = [1 4 2 6 3 8 2]; %% change data here
n = length(data); % count

% Rank (order) for each value: same as rank function in EXCEL
X = data';
for j = 1:length(X(1,:))
[~, order]= sort(X(:,j));
m(order,j) = 1:length(X(:,j));
end
% rank is 'm'
% Exceedance Probability:
for i=1:length(X)
   P(i)=m(i)./(n+1);
end
% Flow Duration Curve Plot
scatter(P,X);
title('Flow Duration Curve');
xlabel('Exceedance Probability');
ylabel('Flow (MCum)');


% The same thing is done in EXCEL
% Reference : http://crk.iri.columbia.edu/water/course/view.php?id=14
