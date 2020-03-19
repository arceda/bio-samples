# A COMPARISON OF KAMERIS AND CASTOR-KRFE

Kameris and Castor-KRFE are methods based on k-mer frequencies for viral subtypying classification. Kameris was proposed by Solis-Reyes et al. (2018), they compute the k-mer frequencies using a Frequency Chaos Game Representation (FCGR) mean while Castor-KRFE,
proposed by Lebatteux, Remita, and Diallo (2019), compute the k-mer frequencies from the whole training datasets.

# Files and folder descripción:
- chaos_game.- Compute the CGR of a viral genome. The implementation is in Javascript.
- hiv1-genomes.- Some genome samples of HIV-1.
- results.- Here, there are CSVs and images computed with the comparison of Kameris and Castor-kRFE.
- feature_extractor.py.- Implementation of Castor-KRFE [1].
- mykameris...py.- Implementation of Kameris [2].


[1] Lebatteux, D., A. M. Remita, and A. B. Diallo. 2019. “Toward an alignment-free method for feature extraction
and accurate classification of viral sequences.” Journal of Computational Biology 26 (6): 519–535.

[2] Solis-Reyes, S., M. Avino, A. Poon, and L. Kari. 2018. “An open-source k-mer based machine learning tool for
fast and accurate subtyping of HIV-1 genomes.” PloS one 13 (11).
