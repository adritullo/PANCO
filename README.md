
#################################################################################
 _______     _       ____  _____   ______    ___                     _   _          
|_   __ \   / \     |_   \|_   _|.' ___  | .'   `.                  (_) / |_        
  | |__) | / _ \      |   \ | | / .'   \_|/  .-.  \  .--.  __   _   __ `| |-'.---.  
  |  ___/ / ___ \     | |\ \| | | |       | |   | | ( (`\][  | | | [  | | | / /__\\ 
 _| |_  _/ /   \ \_  _| |_\   |_\ `.___.'\\  `-'  /  `'.'. | \_/ |, | | | |,| \__., 
|_____||____| |____||_____|\____|`.____ .' `.___.'  [\__) )'.__.'_/[___]\__/ '.__.' 
      
####################################################################################	  

Project PANCO (PANsharpening and COregistration suite) - V.1.0 (06/06/2024)
Adriano Tullo, INAF Astronomical Observatory of Padua (IT)

Every planetary exploration study faces a common challenge: integrating images from diverse sources, each with unique spatial and spectral resolutions, and often, georeferencing inaccuracies. This software suite is dedicated to facilitating the integration of surface images of planetary bodies through simplified co-registration procedures and reciprocal resolution enhancement using CS pansharpening techniques. 
The first validation study is currently being reviewed in Planetary and Space Science, and involves tests on the mutual implementation of ESA TGO CaSSIS and NASA MRO HiRISE images.

Please cite: Tullo, A., Re, C., Cremonese, G., Martellato, E., La Grassa, R., & Thomas, N. (2024). Performance evaluation of pansharpening for planetary exploration: A case study on the implementation of TGO CaSSIS with MRO HiRISE. Planetary and Space Science, 254, 105997. https://doi.org/10.1016/j.pss.2024.105997

Acknowledgements: 
The study has been supported by the Italian Space Agency (ASI-INAF agreement no. 2020-17-HH.0) and INAF (INAF MiniGrant “Combined implementation of CaSSIS and HiRISE data through pansharpening experiments” - CUP C93C23008430001). Thanks to the CaSSIS Team for the precious advice and feedback.

USAGE:
The scripts are designed as modular functions, so they can be used either individually or imported as a function within a script.
For single use, place the input data path in string format at the bottom of the script together with the index of the bands to be used. For CaSSIS cubes, in the case of complete data of all bands, the order is 1:NIR, 2:RED, 3:PAN, 4:BLUE.


