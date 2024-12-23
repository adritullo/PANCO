
#################################################################################
 _______     _       ____  _____   ______    ___                     _   _          
|_   __ \   / \     |_   \|_   _|.' ___  | .'   `.                  (_) / |_        
  | |__) | / _ \      |   \ | | / .'   \_|/  .-.  \  .--.  __   _   __ `| |-'.---.  
  |  ___/ / ___ \     | |\ \| | | |       | |   | | ( (`\][  | | | [  | | | / /__\\ 
 _| |_  _/ /   \ \_  _| |_\   |_\ `.___.'\\  `-'  /  `'.'. | \_/ |, | | | |,| \__., 
|_____||____| |____||_____|\____|`.____ .' `.___.'  [\__) )'.__.'_/[___]\__/ '.__.' 
      
####################################################################################	  

Project PANCO (PANsharpening and COregistration suite)
Adriano Tullo, INAF Osservatorio Astronomico di Padova
V.0.1 - 06/06/2024

Ogni studio di esplorazione planetaria affronta una sfida comune: integrare immagini provenienti da fonti diverse, ciascuna con risoluzioni spaziali e spettrali uniche e spesso con imprecisioni di georeferenziazione. 
Il presente lavoro presenta la prima versione di suite di script dedicata a facilitare l'integrazione delle immagini di superficie dei corpi planetari attraverso procedure di co-registrazione semplificate e il miglioramento della risoluzione reciproca mediante tecniche di CS pansharpening. 
Il primo studio di validazione è attualmente in corso di revisione su Planetary and Space Science, e prevede test sull'implementazione reciproca di immagini ESA TGO CaSSIS e NASA MRO HiRISE.

Referenza: Please cite: Tullo, A., Re, C., Cremonese, G., Martellato, E., La Grassa, R., & Thomas, N. (2024). Performance evaluation of pansharpening for planetary exploration: A case study on the implementation of TGO CaSSIS with MRO HiRISE. Planetary and Space Science, 254, 105997. https://doi.org/10.1016/j.pss.2024.105997


UTILIZZO:
Per motivi tecnici legati al repository gli script sono archiviati in formato testuale, per utilizzarli l'estensione deve essere sostituita con il formato python .py.
Gli script sono progettati come funzioni modulari, per cui possono essere utilizzate sia singolarmente che importate come funzione all'interno di uno script .
Per l'utilizzo singolo porre il path dei dati di input in formato stringa nella parte inferiore dello script insieme all'indice delle bande da utilizzare. Per CaSSIS, nel caso di dato completo di tutte le bande, l'ordine è 1:NIR, 2:RED, 3:PAN, 4:BLU.