%% Notes / To Do
    % > rand() < or rand() <= ?

%% Toggles / Parameters
    % These are the only things that may need to change between runs. 
    % Variables are in camelCase; toggles are in snake_case
  
  % ========================================================== %
    % Plot animated or just end plot?
    plot_anim = false;
  % ========================================================== %
    % Maximum time steps
    tMax = 10000; 
  % ========================================================== %
    % Population count
    n = 1000;
  % ========================================================== %
    % Number of trait loci; possible phenotype range: 0:nLoci + 1
    nLoci = 30; 
  % ========================================================== %
    % Number of preference loci
    nLociP = nLoci; 
  % ========================================================== %
    % Number of bins
    %   If nLoci is divisible by nBins (i.e., bin edges fall on integers),
    %   bins include lower bound and not upper
            % Per MATLAB Docs:
            % "Each bin includes the leading edge, but does not 
            %  include the trailing edge, except for the last bin 
            %  which includes both edges."
    nBins = 5; 
  % ========================================================== %
    % Initial population average for trait locus sum
    initTraitAvg = 0.5;
  % ========================================================== %
    % Initial population average for trait locus sum
    initPrefAvg = 0.5;
  % ========================================================== %
    %Female preference fn sigma 
    sigF =  0.5;
  % ========================================================== %
    % Psychophysical preferences?
    psy_pref_exp = false;  % prefs multiply normal pref dist by e^x
    psy_pref_cdf = false; % prefs change infl. point of cdf
  % ========================================================== %
    % Weak stabilizing selection per Lande 1981
    % Male viability selection fn sigma
    nat_sel = true; 
    sigM = 1;
    
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

%rng(2);


% Troubleshooting param - how many matings are being skipped? 
nSkip = zeros(tMax,1);

% Define bin edges based on nBins
edges = (0:(1/nBins):1);

% Viability fn:
d = makedist('Normal','mu',0.5,'sigma',sigM);
dMax = pdf(d,0.5);
    % Viability at phenotype p (0 < p < 1) is = pdf(d,p) / dMax
    % Normalized by dMax so max is 1 regardless of sigma


% Setup figures
figure
% Define title ahead of time so it is across all subplots
sgtitle(["t = " + tMax + ";" + "   n = " + n + ";" + "   Loci = "+nLoci + ";" + "   Bins = "+ nBins, "Male viability \sigma = "+ sigM + ";" + "  Female preference \sigma = " + sigF])


% Initialize females
% Data structure: 
%   [Trait sum     Pref sum     Trait loci     Preference loci]   
F = [zeros(n,1), zeros(n,1), ...
    (rand(n,nLoci) < initTraitAvg), (rand(n,nLociP) < initPrefAvg)];
% Own trait sum (proportion 0/1) 
F(:,1) = sum(F(:,3:(2+nLoci)),2)  ./ nLoci;
% Preference sum (proportion 0/1) 
F(:,2) = sum(F(:,(3+nLoci):(2+nLoci+nLociP)),2)  ./ nLociP;

% Initialize males
%   [Trait sum     Pref sum     Trait loci     Preference loci]   
M = [zeros(n,1), zeros(n,1), ...
    (rand(n,nLoci) < initTraitAvg), (rand(n,nLociP) < initPrefAvg)];
% Own trait sum (proportion 0/1) 
M(:,1) = sum(M(:,3:(2+nLoci)),2)  ./ nLoci;
% Preference sum (proportion 0/1) 
M(:,2) = sum(M(:,(3+nLoci):(2+nLoci+nLociP)),2)  ./ nLociP;


% Make normally-distributed weights for each bin center to use when females
% choose males from a distribution of bins (not a continuous distribution)
points = (edges(:,1:nBins) + (1/(nBins * 2))).'; % Starts at 1 so we have n-1 points
wts = zeros(nBins); % Initialize array
p = makedist('Normal','mu',0,'sigma',sigF); 


% mu is 0 because we adjust center based on each point
% Calculate PDF(points) for each bin center:
for j = 1:nBins
    if ~psy_pref_cdf
        if ~psy_pref_exp
            wts(:,j) = pdf(p,(points - points(j)));
        else
            wts(:,j) = pdf(p,(points - points(j))) .* exp(points);
        end
    else
        wts(:,j) = cdf(p,(points - points(j)));
    end
end


%% Life cycle

% Timer for benchmarking; only track if not doing animation
if ~plot_anim
    tic
end

% Loop start
for t = 1:tMax
    if nat_sel
        % Delete rows (males) with trait < a unique random value for each male: 
        M(((pdf(d,M(:,1)) ./ dMax) <= rand((size(M,1)),1)),:) = [];
    end
    % Store bin number for males (used in female selection) 
    % (Leading/trailing edge rules of discretize () are as described in intro)
    discM = discretize(M(:,1).', edges).';

    % Mate choice

        % In order to have females choose from bins (categorical) per the 
        % female preference function - i.e., to have her choose around the sum
        % of her own preference loci, following a distribution:
        % e.g. assume 5 bins, the probability that she will choose a male 
        % from a bin is proportional to the bin's height in her fitness function: 
        %
        %      █ 
        %    █ █ █
        % '_'_'_'_'_'
        %  1 2 3 4 5
        %
        % ^ In this case, she is most likely to choose a male from bin 3,
        %   but may instead choose one in 2 or 4; she will not choose 1 or 5. 
        %
        % Female will select from all males with each male weighted 
        % according to the value/weight of its bin in her fitness function.

    % Store bin number for females to match with males, from preference loci
    discF = discretize(F(:,2).', edges ).';
    % Initialize array for storing preference wts for each male,
    % for females in each bin:
    wtsM = zeros(size(M,1), nBins); %need to use size of M because of natsel

    % Not fully preallocated because length is stochastic, but initialized:
    % (stores pairing information)
    O = [];

    for b = 1:nBins
      if ~sum(discF==b) %if no females in bin, move to next bin
        nSkip(t,1) = nSkip(t,1) +1; % For troubleshooting - store that matings were skipped
        continue
      end  

      % From discM, need to calc weights from  for each male index, per bin
      % The following replaces the discretized "bin numbers" in discM,
      % preserving index, with the associated weights for that bin given
      % the current (female) bin:
      wtsM = changem(discM,wts(:,b),1:nBins); 
      % These weights can be passed to datasample so she selects from all
      % males, each with an associated weight based on the preference
      % function and her preferred bin.
            
      % Using datasample(), pull random rows from males, using weights;
      % WITH REPLACEMENT (polygyny)
      % with output array size rows = (nFemales with bin b)
      % Add male rows to the rows for the females in question -
      % see explanation of haploid mating in Of, Om below.
      % Add the new rows to the bottom of the old rows prom previous bins (O;)
      O = [O; datasample(M, sum(discF==b), 'Weights', wtsM)...
        + F(discF==b,:)];
      % O stores the PAIRINGS, from which the OFFSPRING will be randomly
      % extracted below in Of, Om.
    end

    % Offspring (O) may not be as long as F or M  if no males in bin
    % Draw from pairings (O) randomly until count reaches n; males and females
        Of = O(randsample(1:size(O,1),n,true).',:);
        Om = O(randsample(1:size(O,1),n,true).',:);

    % HAPLOID MATING
    % This is done in a slightly unusual way for performance & to avoid loops.
    % It assumes free recombination (r = 0.5).
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
    Of(:,3:nLoci+nLociP+2) = (rand(size(Of(:,3:nLoci+nLociP+2))) <= (Of(:,3:nLoci+nLociP+2)./2));
    % Store sum of own trait loci: (proportion 0/1) 
    Of(:,1) = sum(Of(:,3:(2+nLoci)),2)  ./ nLoci;
    % Store sum of own preference loci: (proportion 0/1) 
    Of(:,2) = sum(Of(:,(3+nLoci):(2+nLoci+nLociP)),2)  ./ nLociP;

    % Repeat for male offspring: 
    Om(:,3:nLoci+nLociP+2) = (rand(size(Om(:,3:nLoci+nLociP+2))) <= (Om(:,3:nLoci+nLociP+2)./2));
    %Store sum of own trait loci: (proportion 0/1) 
    Om(:,1) = sum(Om(:,3:(2+nLoci)),2)  ./ nLoci;
    % Store sum of own preference loci: (proportion 0/1) 
    Om(:,2) = sum(Om(:,(3+nLoci):(2+nLoci+nLociP)),2)  ./ nLociP;
    
    % Pass to F/M arrays for next generation:
    F = Of;
    M = Om;

    % Timer for benchmarking; only track if not doing animation
    if ~plot_anim && t == tMax
        toc
    end

    % If animated OR if last timestep - show plots
    if plot_anim | t == tMax 
        pause(0.1) % Animation speed
        
        % Setup
        % How MATLAB handles values at bin edges: 
        % "Each bin includes the leading edge, but does not 
        %  include the trailing edge, except for the last bin 
        %  which includes both edges."

        % Need to define edges of histogram 'bins' for locus count
        plotbins = ((0:1:nLoci))./nLoci;

        % For plotting generic female preference fn shape: define per sig_f
        p = makedist('Normal','mu',0.5,'sigma',sigF);
        xs = 0:0.05:1;
        ys = pdf(p,0:0.05:1);
        pMax = pdf(p,0.5);

        % if ~psy_pref_cdf
        %     if ~psy_pref_exp
        % 
        %     else
        %         ys = pdf(p,0:0.05:1) .* exp(xs);
        %         pMax = max(ys);
        %     end
        % 
        % else
        %     ys = cdf(p,0:0.05:1);
        %     pMax = cdf(p,1);
        % end
        
        % Actual plotting:
        % Males
        subplot(1,2,1)
        % Plot viability dist
        plot(0:0.05:1,pdf(d,0:0.05:1)./dMax)
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
        plot(xs,ys./pMax, Color='#eb9866', LineStyle=':')
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
