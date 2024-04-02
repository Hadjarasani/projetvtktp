# projetvtktp
EXPLICATION DU CODE :
J’ai repris les techniques de l’out-of-core du tp7 avec le parallélisme vu en cours pour optimiser mon code au maximum et avoir un bon rendu des images. Globalement, voici comment fonctionne mon code:

Dans un premier temps je répartis les Z entre chaque processus;
Chaque processus divise ses Z par numPasses et charge une partie de ses Z dans le fichier raw;
Ensuite je fais une boucle pour effectuer plusieurs azimuth et ou modifier modifier les valeurs de coupe à chercher;
Chaque processus fait le rendu de ses Z en plusieurs fois dans le sens où il charge les données pour une partie de ses Z, fait le rendu et garde le rendu final pour ses Z via un zbuffer local;
Enfin je fait un MPI_Gather sur root pour récupérer tous les rgba et zbuffer et calculer l’image finale, une fois que tous les processus ont fait le rendu de leurs Z. Pour avoir des images de meilleure qualité, j’ai augmenté la taille de la fenêtre de rendu à 1400.

Cette approche permet de charger les images en exploitant pleinement les capacités de chaque cœur de processeur, tout en minimisant les besoins en mémoire vive, même pour des fichiers bruts de grande taille. Cependant, un inconvénient majeur réside dans le fait que chaque processus MPI effectue plusieurs lectures du fichier brut pour générer ses propres données RGBA locales, avant de les fusionner pour obtenir l'image RGBA finale.
