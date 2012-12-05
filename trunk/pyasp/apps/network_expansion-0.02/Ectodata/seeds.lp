%seed("_2__45__KETOGLUTARATE").						% A ENLEVER. INTERMEDIAIRE. ECTOCARPUS est sense le produire.
seed("Acceptor"). 							% OK. Accepteur d'electrons oxydé. Générique qui intervient dans beaucoup de réactions.
%seed("ACETYL__45__COA").  						% A ENLEVER. INTERMEDIAIRE DE REACTION. 
seed("ACP").								% OK ? Proteine produite pas plusieurs genes chez ectocarpus. Quel est son role dans le réseau ? Nécessaire à des réactions métaboliques et du transport. On doit supposer qu'elle est présente et préalablement produite par la machinerie cellulaire. Vérifier si elle utilisée.
%seed("ADENOSYL__45__HOMO__45__CYS").					% A ENLEVER. AA produit par ectocarpus.
seed("ADP").								% OK.	
seed("AMMONIA").							% OK.
seed("AMP").								% OK	
seed("ATP").								% OK	
seed("BIOTIN").								% OK
seed("CA__43__2").		% calcium				% OK 	
seed("CARBON__45__DIOXIDE").						% OK
seed("CDP__45__12921"). %nitrate de sodium  				% OK
seed("CL__45__").							% OK
seed("CO__45__A").							% OK
seed("CPD0__45__1695"). %chlorure de cobalt   				% OK	
seed("CPD0__45__2394"). %chlorure de manganese				% OK
seed("CPD__45__3"). 	%molybdate					% OK
seed("CPD__45__315"). %cyanocob(III)alamin				% OK
seed("CU__43__").							% OK
seed("CU__43__2").							% OK	
seed("Donor__45__H2").							% OK. Ca va avec "Accceptor"
seed("EDTA").								% OK
seed("GDP"). 								% OK
%seed("GLN").								% A ENLEVER. Synthetise a partir de l'ammoniac.		
%seed("GLT").								% A ENLEVER. Synthetise a partir de l'ammoniac. ? Montrer que un des deux est nécessaire à cause d'un effet de cycles.		
%seed("GLUTATHIONE").  							% A ENLEVER. INTERMEDIAIRE Produits par ectocarpus.
seed("GTP").								% OK	
seed("HYDROGEN__45__PEROXIDE"). 					% OK COFACTEUR
seed("K__43__").							% OK						
%seed("MALONYL__45__ACP").						% A ENLEVER. Produit par ectocarpus. Risque de soucis à cause des cycles.
seed("MG__43__2").							% OK
seed("NA__43__").							% OK	
seed("NAD__45__P__45__OR__45__NOP").					% OK
seed("NAD").								% OK
seed("NADH__45__P__45__OR__45__NOP").					% OK
seed("NADH").								% OK					
seed("NADP").								% OK
seed("NADPH").								% OK	
seed("NITRATE").							% OK
seed("OXYGEN__45__MOLECULE").						% OK
seed("Pi").								% OK
seed("PPI").								% OK
seed("PROTON").								% OK
%seed("PYRUVATE").							% A ENLEVER. Risque de soucis à cause des cycles. IL doit falloir un composé carboné pour initialiser le système
%seed("S__45__ADENOSYLMETHIONINE").					% A ENLEVER 	
seed("SULFATE").							% OK
seed("THIAMINE").							% OK
seed("UDP").								% OK  Pourquoi il ne faut pas de l'UTP ? Sans doute symétriser le tout en rajoutant l'UTP
seed("Vitamins__45__B12").						% OK
seed("COB__45__I__45__ALAMIN").						% OK vit B12	
seed("WATER").								% OK
seed("ZN__43__2").							% OK
