%% Notes / To Do
% > 


%% Toggles / Parameters
% These are the only things that may need to change between runs. 

% Plot animated or just end plot?
plotAnim = false;

% Maximum time steps
tmax = 1000; 

% Population count
n = 1000; 

% Number of trait loci; total number of possible phen. = nLoci + 1 (includes 0)
nLoci = 10; 

% Number of preference loci
nLoci_P = nLoci; 

% Number of bins
%   If nLoci is divisible by nBins (i.e., bin edges fall on integers),
%   bins include lower bound and not upper
        % Per MATLAB Docs: 
        % "Each bin includes the leading edge, but does not 
        %  include the trailing edge, except for the last bin 
        %  which includes both edges."
nBins = 5; 

%Female preference fn sigma 
sig_f =  0.1;

% Weak stabilizing selection per Lande 1981
% Male viability selection fn sigma
natsel = true; 
sig_m = 1;

    % To view fitness fn shape: 
    
    % plot(d) %not normalized4figure
    % figure
    % plot([0:0.05:1],pdf(d,[0:0.05:1])./dMax) %plot fitness dist
    % xlim([0,1])
    % ylim([0,1.05])
    % title("Fitness function")


%% Other setup
% RNG seed - un-commenting & setting this will give you reproducible results, i.e. results
% with the same code using the same seed should be identical despite stochasticity
% From MATLAB docs: 
    %Specify seed as a nonnegative integer, such as rng(1), to initialize the random number generator with that seed.
    %Specify seed as "shuffle" to initialize the generator seed based on the current time.

%rng(1);


% Troubleshooting param - how many matings are being skipped? 
nSkip = zeros(tmax,1);

% Define bin edges based on nBins
edges = (0:(1/nBins):1);

% Viability fn:
d = makedist('Normal','mu',0.5,'sigma',sig_m);
dMax = pdf(d,0.5);
    % Viability at phenotype p (0 < p < 1) is = pdf(d,p) / dMax
    % Normalized by dMax so max is 1 regardless of sigma


% Setup figures
figure
% Define title ahead of time so it is across all subplots
sgtitle(["t = " + tmax + ";" + "   n = " + n + ";" + "   Loci = "+nLoci + ";" + "   Bins = "+ nBins, "Male viability \sigma = "+ sig_m + ";" + "  Female preference \sigma = " + sig_f])

% Initialize females
% Data structure: 
%   [Pref. sum     Own sum     Preference loci     Trait loci]   
F = [zeros(n,1), zeros(n,1), randi([0,1],[n,nLoci + nLoci_P])];
% Own trait sum (proportion 0/1) 
F(:,1) = sum(F(:,3:(2+nLoci)),2)  ./ nLoci;
% Preference sum (proportion 0/1) 
F(:,2) = sum(F(:,(3+nLoci):(2+nLoci+nLoci_P)),2)  ./ nLoci_P;

% Initialize males
%   [Pref. sum     Own sum     Preference loci     Trait loci] 
M = [zeros(n,1), zeros(n,1), randi([0,1],[n,nLoci + nLoci_P])];
% Own trait sum (proportion 0/1) 
M(:,1) = sum(M(:,3:(2+nLoci)),2)  ./ nLoci;
% Preference sum (proportion 0/1) 
M(:,2) = sum(M(:,(3+nLoci):(2+nLoci+nLoci_P)),2)  ./ nLoci_P;


% Make normally-distributed weights for each bin center to use when females
% choose males from a distribution of bins (not a continuous distribution)
points = (edges(:,1:nBins) + (1/(nBins * 2))).'; % Starts at 1 so we have n-1 points
wts = zeros(nBins); % Initialize array
p = makedist('Normal','mu',0,'sigma',sig_f); 
% mu is 0 because we adjust center based on each point
% Calculate PDF(points) for each bin center:
for j = 1:nBins
    wts(:,j) = pdf(p,(points - points(j)));
end


%% Life cycle

% Timer for benchmarking; only track if not doing animation
if ~plotAnim
    tic
end

% Loop start
for t = 1:tmax
    if natsel
        % Delete rows (males) with trait < a unique random value for each male: 
        M(((pdf(d,M(:,1)) ./ dMax) <= rand((size(M,1)),1)),:) = [];
    end
    % Store bin number for males (used in female selection) 
    % (Leading/trailing edge rules of discretize () are as described in intro)
    discM = discretize(M(:,1).', edges).';

    % Mate choice

    %   Assume polygyny / replacement

        % In order to have females choose from bins (categorical) per the 
        % female preference function - i.e., to have her choose around the sum
        % of her own preference loci, following a distribution:
        % e.g. assume 5 bins, the probability that she will choose from a bin 
        % is proportional to its height in her fitness function: 
        %
        %      █ 
        %    █ █ █
        % '_'_'_'_'_'
        %  1 2 3 4 5
        %
        % ^ In this case, she is most likely to choose to search within bin 3,
        %   but may instead search in 2 or 4; she will not search in 1 or 5. 
        %
        % The way this was implemented was to assign her a bin according to 
        % the distribution, then choose randomly in that bin. 
        % This is INSTEAD of having her select from all males with each male
        % weighted according to the value of its own bin in her fitness
        % function, because - in short - doing that with datasample() 
        % produced different distributions less clearly tied to sig_m.  

    % Store bin number for females to match with males
    discF = discretize(F(:,2).', edges ).';
    % Assign females bins based on mating weight fns
    % The following could also just be [1:n]; it uses discretize for
    % consistency (especially with leading/tailing edge problems)
    points_d = discretize(points.', edges ).'; 
    binsF = zeros(n,1);

    % Female preference fn adjusts bin of choice per sig_f distribution 
    for b = 1:nBins
      if ~sum(discF==b) %if no females in bin, move to next bin
        nSkip(t,1) = nSkip(t,1) +1; % For troubleshooting - store that matings were skipped
        continue
      end  
      % The bin that a female chooses is adjusted FROM her default bin (at
      % the mean of her distribution) TO an adjacent bin, weighted based on
      % the preference function evaluated at each bin location. 
      % Most likely to stay at center bin - might be adjusted up/down. 
      binsF(discF==b,1) = datasample(points_d,sum(discF==b),'Weights', wts(:,b).');
      % ^ This returns (n = how many females have the bin in question)
      % randomly selected values from the full list of 'points' of each bin,
      % weighted by the fitness function that would be defined using the
      % bin in question as the center. 
    end
    

    %Not fully preallocated because length is stochastic, but initialized:
    O = [];
    
    % Index mask of females in bin b = (binsF==b)
    for b = 1:nBins
        % Skip if no females in this bin:
        if ~sum(binsF==b) || ~sum(discM==b) 
        continue
        end
    
    % Using datasample(), pull random row (index) from group "males of type b"
    % WITH REPLACEMENT (polygyny)
    % with array size 1x(nFemales with bin b)
    % Call that index from M to get rows for each male
    % Add that male row to the rows for the females in question -
    % see explanation of haploid mating in Of, Om
    O = [O; M(datasample(find(discM == b),sum(binsF==b)),:) + F(binsF==b,:)];
    end

    % Offspring (O) may not be as long as F or M  if no males in bin
    % Draw from pairings randomly until count reaches n; males and females
        Of = O(randsample(1:size(O,1),n,true).',:);
        Om = O(randsample(1:size(O,1),n,true).',:);

    % HAPLOID MATING
    % This is done in a slightly unusual way for performance & avoiding loops.
    % It does assume free recombination (r = 0.5).
    % The binary loci from the male and female are added together, creating
    % values of 0, 1, or 2. 
    % Divide the values by 2; now values are:
        % For mating 0x0 at locus = 0
        % For mating 1x0 at locus = 0.5
        % For mating 1x1 at locus = 1
    % These values are used as the probability that the offspring's locus 
    % will be filled with a 1 instead of 0. 
    % 1x1 matings are always filled with 1; 0x0 always with 0; 1x0 are
    % decided per a 50% chance. 

    % Female offspring
    % Assign these loci as above: 
    Of(:,3:nLoci+nLoci_P+2) = (rand(size(Of(:,3:nLoci+nLoci_P+2))) <= (Of(:,3:nLoci+nLoci_P+2)./2));
    % Store sum of own trait loci: (proportion 0/1) 
    Of(:,1) = sum(Of(:,3:(2+nLoci)),2)  ./ nLoci;
    % Store sum of own preference loci: (proportion 0/1) 
    Of(:,2) = sum(Of(:,(3+nLoci):(2+nLoci+nLoci_P)),2)  ./ nLoci_P;

    % Repeat for male offspring: 
    Om(:,3:nLoci+nLoci_P+2) = (rand(size(Om(:,3:nLoci+nLoci_P+2))) <= (Om(:,3:nLoci+nLoci_P+2)./2));
    %Store sum of own trait loci: (proportion 0/1) 
    Om(:,1) = sum(Om(:,3:(2+nLoci)),2)  ./ nLoci;
    % Store sum of own preference loci: (proportion 0/1) 
    Om(:,2) = sum(Om(:,(3+nLoci):(2+nLoci+nLoci_P)),2)  ./ nLoci_P;
    
    % Pass to F/M arrays for next generation:
    F = Of;
    M = Om;

    % Timer for benchmarking; only track if not doing animation
    if ~plotAnim && t == tmax
        toc
    end

    % If animated OR if last timestep - show plots
    if plotAnim | t == tmax 
        pause(0.1) % Animation speed
        
        % Setup

        % How MATLAB handles values at bin edges: 
        % "Each bin includes the leading edge, but does not 
        %  include the trailing edge, except for the last bin 
        %  which includes both edges."

        % Need to define edges of histogram 'bins' for locus count
        plotbins = ((0:1:nLoci))./nLoci;

        % For plotting generic female preference fn shape: define per sig_f
        p = makedist('Normal','mu',0.5,'sigma',sig_f);
        pMax = pdf(p,0.5);

        % Actual plotting:
        % Males
        subplot(1,2,1)
        % Plot viability dist
        plot([0:0.05:1],pdf(d,[0:0.05:1])./dMax)
        hold on
        % Plot bins - where each male falls in the population-defined bins
        % * different from female plot
        histogram(M(:,1), edges, Normalization='probability',...
            EdgeColor='#1171be', FaceColor='none')
        % Plot alleles - sum of own trait loci
        histogram(M(:,1), plotbins, Normalization='probability',...
            EdgeColor='k', FaceColor='#1171be')
        xlim([0,1])
        ylim([0,1.05])
        title("Male")
        hold off
    
        % Females
        subplot(1,2,2)
        % Plot generic preference fn shape for females
        plot([0:0.05:1],pdf(p,[0:0.05:1])./pMax, Color='#eb9866', LineStyle=':')
        hold on
        % Plot bins - !different from males! - sum of each female's
        % preference loci. Different color helps remember
        histogram(F(:,2), edges, Normalization='probability',...
            EdgeColor='k', FaceColor='none')
        % Plot alleles - sum of own trait loci
        histogram(F(:,1), plotbins, Normalization='probability',...
            EdgeColor='k', FaceColor='#dd5400')
        xlim([0,1])
        ylim([0,1.05])
        title("Female")
        hold off
    end
end


