% How to call
% gringo --const n=1 promoter_facts.lp correlationelation_data.lp consistency.lp | clasp 0 --opt-all --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlationelation_matrix.lp consistency.lp | clasp 0 --opt-heu --project
% gringo --const n=3 oxi_tfbs.lp oxi_correlationelation_matrix.lp consistency.lp | clasp 0 --opt-heu --project --restart-on-model --heu=vsids --restarts=2
%gringo consistency.lp --const n=3 ../../../test/oxi_promo/oxi_tfbs.lp ../../../test/oxi_promo/oxi_correlation_matrix.lp | /home/wv/bin/linux/32/clasp-mt-banane 0 --opt-heu --project --restart-on-model --heu=vsids --restarts=2
% when have two anti-/corellated genes a tf incommon
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), correlation(GENA,GENB,1).
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), correlation(GENA,GENB,-1).

% when have have two uncorrelationelated genes a distinct tf
one_different_tf(GENA,GENB) :- tf(GENA,TF), not tf(GENB,TF), correlation(GENA,GENB,0).
one_different_tf(GENA,GENB) :- tf(GENB,TF), not tf(GENA,TF), correlation(GENA,GENB,0).

% when have two uncorrelationelated genes all tf's incommon
all_common_tf(GENA,GENB) :- not one_different_tf(GENA,GENB), correlation(GENA,GENB,0).


%%%% define a repair

% number of added unknown transcription factors
unknowntf(1..n).

% define possible repairs on regulation model
%pos(add_tf_both(GENA,GENB,TF)) :- correlation(GENA,GENB, 1), not one_common_tf(GENA,GENB), unknowntf(TF).
%pos(add_tf_both(GENA,GENB,TF)) :- correlation(GENA,GENB,-1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENA,TF)) :- correlation(GENA,GENB, 1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENB,TF)) :- correlation(GENA,GENB, 1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENA,TF)) :- correlation(GENA,GENB,-1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENB,TF)) :- correlation(GENA,GENB,-1), not one_common_tf(GENA,GENB), unknowntf(TF).
pos(add_tf(GENA,TF)) :- correlation(GENA,GENB,0),      all_common_tf(GENA,GENB), unknowntf(TF) .
pos(add_tf(GENB,TF)) :- correlation(GENA,GENB,0),      all_common_tf(GENA,GENB), unknowntf(TF) .
pos(rm(TF)) :- tf(GENA,TF).

% choose applied repairs
{app(REP)} :- pos(REP).

% minimize applied repairs
#minimize [ not false=0, app(R) : pos(R) ].

% effects of repair
tf_ur(GENE,TF) :- app(add_tf(GENE,TF)).
%tf_ur(GENA,TF) :- app(add_tf_both(GENA,GENB,TF)).
%tf_ur(GENB,TF) :- app(add_tf_both(GENA,GENB,TF)).
tf_ur(GENE,TF) :- tf(GENE,TF), not app(rm(TF)).
%tf_ur(GENE,TF) :- tf(GENE,TF).

%%%% define model under minimal repair

% when have under repair two genes a tf in common
one_common_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), tf_ur(GENB,TF), correlation(GENA,GENB,1).
one_common_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), tf_ur(GENB,TF), correlation(GENA,GENB,-1).

% when have under repair two genes a distinct tf
one_different_tf_ur(GENA,GENB) :- tf_ur(GENA,TF), not tf_ur(GENB,TF), correlation(GENA,GENB,0).
one_different_tf_ur(GENA,GENB) :- tf_ur(GENB,TF), not tf_ur(GENA,TF), correlation(GENA,GENB,0).

% when have under repair two genes all tf's incommon
all_common_tf_ur(GENA,GENB) :- not one_different_tf_ur(GENA,GENB), correlation(GENA,GENB,0).


%%%% define consistency of {repaired} model


% when are correlationelation data and gene regulation model inconsistent
:- correlation(GENA,GENB, 1), not one_common_tf_ur(GENA,GENB).
:- correlation(GENA,GENB, 0),     all_common_tf_ur(GENA,GENB).
:- correlation(GENA,GENB,-1), not one_common_tf_ur(GENA,GENB).
 

#hide correlation/3.
#hide tf/2.
#hide tf_ur/2.

#hide.
#show app/1.
%#show conflict/2.

gene(X)  :- correlation(X,Y,Z).
gene(Y)  :- correlation(X,Y,Z).

gnFirst(X) :- gene(X), gene(Y), Y < X.
gFirst(X)  :- not gnFirst(X), gene(X).

gnLast(X) :- gene(X), gene(Y), Y > X.
gLast(X)  :- not gnLast(X), gene(X).

gnNext(X,Z) :- gene(X;Y;Z), X < Y, Y < Z, X < Z.
gNext(X,Y)  :- gene(X;Y), X < Y, not gnNext(X,Y).

gene(X,TF)  :- tf_ur(X,TF), unknowntf(TF).
gene(Y,TF)  :- gNext(X,Y), gene(X,TF), unknowntf(TF).
:- gene(X,TF+1), not gene(X,TF), unknowntf(TF;TF+1).









