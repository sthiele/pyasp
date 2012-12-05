% How to call
% gringo --const n=1 promoter_facts.lp correlation_data.lp consistency.lp | clasp 0 --opt-all --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlation_matrix.lp consistency.lp | clasp 0 --opt-heu --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlation_matrix.lp consistency.lp | clasp 0 --opt-heu --project --restart-on-model --heu=vsids --restarts=2

%%% compute the full correlation matrix
% guess a correlation for unsure (q1, q2) measurements
% 1{gcorrelation(GENA,GENB, 1), gcorrelation(GENA,GENB,0)}1 :- correlation(GENA,GENB,q1).
% 1{gcorrelation(GENA,GENB,-1), gcorrelation(GENA,GENB,0)}1 :- correlation(GENA,GENB,q2).
% :- gcorrelation(GENA,GENB,X), gcorrelation(GENB,GENA,Y), X!=Y.

% if measurements are sure (1,0,-1) use them
gcorrelation(GENA,GENB,-1) :- correlation(GENA,GENB,-1).
gcorrelation(GENA,GENB, 0) :- correlation(GENA,GENB, 0).
gcorrelation(GENA,GENB, 1) :- correlation(GENA,GENB, 1).


% when have two anti-/corellated genes a tf incommon
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), gcorrelation(GENA,GENB,1).
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), gcorrelation(GENA,GENB,-1).

% when have have two uncorrelated genes a distinct tf
one_different_tf(GENA,GENB) :- tf(GENA,TF), not tf(GENB,TF), gcorrelation(GENA,GENB,0).
one_different_tf(GENA,GENB) :- tf(GENB,TF), not tf(GENA,TF), gcorrelation(GENA,GENB,0).

% when have two uncorrelated genes all tf's incommon
all_common_tf(GENA,GENB) :- not one_different_tf(GENA,GENB), gcorrelation(GENA,GENB,0).


%%%% define a repair

% number of added unknown transcription factors
unknowntf(1..n).

% define possible repairs on regulation model
pos(add_tf_both(GENA,GENB,TF)) :- gcorrelation(GENA,GENB, 1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf_both(GENA,GENB,TF)) :- gcorrelation(GENA,GENB,-1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENA,TF)) :- all_common_tf(GENA,GENB), gcorrelation(GENA,GENB,0), unknowntf(TF) .
pos(add_tf(GENB,TF)) :- all_common_tf(GENA,GENB), gcorrelation(GENA,GENB,0), unknowntf(TF) .
% choose applied repairs
{app(REP)} :- pos(REP).

% minimize applied repairs
#minimize [ app(R) : pos(R) ].

% effects of repair
tf_ur(GENE,TF) :- app(add_tf(GENE,TF)).
tf_ur(GENA,TF) :- app(add_tf_both(GENA,GENB,TF)).
tf_ur(GENB,TF) :- app(add_tf_both(GENA,GENB,TF)).
tf_ur(GENE,TF) :- tf(GENE,TF).


%%%% define model under minimal repair

% when have under repair two genes a tf in common
one_common_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), tf_ur(GENB,TF), gcorrelation(GENA,GENB,1).
one_common_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), tf_ur(GENB,TF), gcorrelation(GENA,GENB,-1).

% when have under repair two genes a distinct tf
one_different_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), not tf_ur(GENB,TF), gcorrelation(GENA,GENB,0).
one_different_tf_ur(GENA,GENB) :- tf_ur(GENB,TF), not tf_ur(GENA,TF), gcorrelation(GENA,GENB,0).

% when have under repair two genes all tf's incommon
all_common_tf_ur(GENA,GENB) :- not one_different_tf_ur(GENA,GENB), gcorrelation(GENA,GENB,0).


%%%% define consistency of {repaired} model


% when are correlation data and gene regulation model inconsistent
%inconsistent(GENA,GENB) :- gcorrelation(GENA,GENB,1), not both(GENA,GENB).
%inconsistent(GENA,GENB) :- gcorrelation(GENA,GENB,0), allsame(GENA,GENB).
 :- gcorrelation(GENA,GENB,1), not one_common_tf_ur(GENA,GENB).
 :- gcorrelation(GENA,GENB,0), all_common_tf_ur(GENA,GENB).
 :- gcorrelation(GENA,GENB,-1), not one_common_tf_ur(GENA,GENB).
 

%#hide pos/1.
#hide gcorrelation/3.
#hide tf/2.
#hide tf_ur/2.

 #hide.
%#show one_common_tf/2.
%#show one_different_tf/2.
#show app/1.