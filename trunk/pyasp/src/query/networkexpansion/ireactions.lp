% inputs from the met_net
% species(Species).
% reaction(Reaction,Network).
% reactant(Species,Reaction).
% product(Species,Reaction).
% seed(Species).
% target(Species).

% find minimal set of xreactions such that xscope contains targets in pscope

%draft(draftnetwork).

% what is producible by the original network

  dscope(M) :- seed(M).

  dscope(M) :- product(M,R,N), reaction(R,N), draft(N),
               dscope(M2) : reactant(M2,R,N).


% what is producible with all reactions

  pscope(M) :- seed(M).

  pscope(M) :- product(M,R,N), reaction(R,N),
               pscope(M2) : reactant(M2,R,N).

% what are the interesting reactions

  interesting(M) :- target(M), not dscope(M), pscope(M).

  interesting(M) :- reactant(M,R,N), iireaction(R,N), not dscope(M).

  iireaction(R,N) :- interesting(M), product(M,R,N), reaction(R,N).


  binteresting(M) :- dscope(M).

  binteresting(M) :- product(M,R,N), bireaction(R,N).

  bireaction(R,N) :- binteresting(M), reactant(M,R,N), product(M2,R,N), reaction(R,N), not dscope(M2).

ireaction(R,N) :- iireaction(R,N), bireaction(R,N), not draft(N).



#hide.
#show ireaction/2.

