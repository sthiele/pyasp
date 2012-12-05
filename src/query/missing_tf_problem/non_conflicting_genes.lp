% How to call
% gringo --const n=1 promoter_facts.lp correlationelation_data.lp consistency.lp | clasp 0 --opt-all --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlationelation_matrix.lp consistency.lp | clasp 0 --opt-heu --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlationelation_matrix.lp consistency.lp | clasp 0 --opt-heu --project --restart-on-model --heu=vsids --restarts=2
%gringo consistency.lp --const n=3 ../../../test/missing_tf_problem/oxi_tfbs.lp ../../../test/missing_tf_problem/oxi_correlation_matrix.lp | /home/wv/bin/linux/32/clasp-mt-banane 0 -t 16 --opt-heu --project --solu --reset-restarts --local-restarts --heu=vsids --restarts=16


% when have two anti-/corellated genes a tf incommon
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), correlation(GENA,GENB,1).
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), correlation(GENA,GENB,-1).

% when have have two uncorrelationelated genes a distinct tf
one_different_tf(GENA,GENB) :- tf(GENA,TF), not tf(GENB,TF), correlation(GENA,GENB,0).
one_different_tf(GENA,GENB) :- tf(GENB,TF), not tf(GENA,TF), correlation(GENA,GENB,0).

% when have two uncorrelationelated genes all tf's incommon
all_common_tf(GENA,GENB) :- not one_different_tf(GENA,GENB), correlation(GENA,GENB,0).


%%%% define conflicts of model


% when are correlationelation data and gene regulation model inconsistent
conflict(GENA) :- correlation(GENA,GENB, 1), not one_common_tf(GENA,GENB).
conflict(GENB) :- correlation(GENA,GENB, 1), not one_common_tf(GENA,GENB).

conflict(GENA) :- correlation(GENA,GENB, 0),     all_common_tf(GENA,GENB).
conflict(GENB) :- correlation(GENA,GENB, 0),     all_common_tf(GENA,GENB).

conflict(GENA) :- correlation(GENA,GENB,-1), not one_common_tf(GENA,GENB).
conflict(GENB) :- correlation(GENA,GENB,-1), not one_common_tf(GENA,GENB).


no_conflict(GENA) :- correlation(GENA,GENB,S), not conflict(GENA).
no_conflict(GENB) :- correlation(GENA,GENB,S), not conflict(GENB).

#hide.
#show no_conflict/1.








