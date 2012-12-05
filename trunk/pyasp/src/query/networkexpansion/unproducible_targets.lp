% inputs from the met_net
% species(M).
% reaction(R).
% reactant(M,R).
% product(M,R).

% seed(M).
% target(M).


% what is producible by all reactions

  pscope(M) :- seed(M).

  pscope(M) :- product(M,R,N), reaction(R,N),
               pscope(M2) : reactant(M2,R,N).

%version roland
%pscope(M) :- seed(M).
%
%_reactant(#infimum,R)  :- reaction(R,_).
%_reactant(#supremum,R) :- reaction(R,_).
%_reactant(M,R) :- reactant(M,R).

%order_reactant(M1,M3,R) :- 
%	_reactant(M1,R), _reactant(M3,R), M1 < M3, 
%	not _reactant(M2,R) : _reactant(M2,R) : M1 < M2 : M2 < M3.

%pscope(M) :- product(M,R), reaction(R,N), has_pscope(M2,R), order_reactant(M2,#supremum,R).
%has_pscope(#infimum,R)  :- _reactant(#infimum,R).
%has_pscope(M2,R)        :- has_pscope(M1,R), order_reactant(M1,M2,R), pscope(M2).



% show interesting predicates

  unproducible_target(M) :- target(M), not pscope(M).



#hide.
#show unproducible_target/1.
