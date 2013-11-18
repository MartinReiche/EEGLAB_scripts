% Perform one way repeated measures analysis of variance RMANOVA. Takes
% Data in matrix form where lines represent the subjects with the data of
% severeal factors represented by columns.
% 
% The scripts are based on - RMAOV1 - Trujillo-Ortiz, A., R. Hernandez-Walls and
% R.A. Trujillo-Perez. (2004). RMAOV1:One-way repeated measures ANOVA. A
% MATLAB file. (http://www.mathworks.com/matlabcentral/fileexchange)
%
% Copyright (c) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function P1 = OneWayrmAoV(testData,alpha)
% evaluate input arguments and set defaults

if nargin < 2,
   alpha = 0.05; %(default)
end; 

if (alpha <= 0 | alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

if nargin < 1, 
   error('Requires at least one input argument.');
   return;
end;

% change the form of the test data to resemble R-like data.frame structure
% with the dependent variable in column 1, the independant variable in
% column 2 and the subject index in column 3

factors = repmat(find(testData(1,:)),size(testData,1),1);
X = [testData(:) factors(:) repmat(1:size(testData,1),1,size(testData,2))'];

k = max(X(:,2));
s = max(X(:,3));


%Analysis of Variance Procedure.
m=[];n=[];nn=[];A=[];
indice = X(:,2);
for i = 1:k
   Xe = find(indice==i);
   eval(['X' num2str(i) '=X(Xe,1);']);
   eval(['m' num2str(i) '=mean(X' num2str(i) ');'])
   eval(['n' num2str(i) '=length(X' num2str(i) ') ;'])
   eval(['nn' num2str(i) '=(length(X' num2str(i) ').^2);'])
   eval(['xm = m' num2str(i) ';'])
   eval(['xn = n' num2str(i) ';'])
   eval(['xnn = nn' num2str(i) ';'])
   eval(['x =(sum(X' num2str(i) ').^2)/(n' num2str(i) ');']);
   m=[m;xm];n=[n;xn];nn=[nn,xnn];A=[A,x];
end;

S=[];
indice=X(:,3);
for j=1:s
   Xe=find(indice==j);
   eval(['S' num2str(j) '=X(Xe,1);']);
   eval(['x =((sum(S' num2str(j) ').^2)/length(S' num2str(j) '));']);
   S=[S,x]; 
end;

C = (sum(X(:,1)))^2/length(X(:,1)); %correction term
SST = sum(X(:,1).^2)-C; %total sum of squares
dfT = length(X(:,1))-1; %total degrees of freedom

SSA = sum(A)-C; %IV sum of squares
v1 = k-1; %IV degrees of freedom
SSS = sum(S)-C; %within-subjects sum of squares
v2 = s-1; %within-subjects degrees of freedom
SSE = SST-SSA-SSS; %error sum of squares
v3 = v1*v2; %error degrees of freedom
MSA = SSA/v1; %IV mean squares
MSS = SSS/v2; %within-subjects mean squares
MSE = SSE/v3; %error mean squares
F1 = MSA/MSE; %IV F-statistic
F2 = MSS/MSE; %within-subjects F-statistic

%Probability associated to the F-statistics.
P1 = 1 - fcdf(F1,v1,v3); 
P2 = 1 - fcdf(F2,v2,v3);   

eta2 = SSA/(SSA+SSE)*100;





