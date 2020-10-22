# network-drawing

Script with different ways to draw monolayer and multilayer networks in R.

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Marco Mello & Renata Muylaert.

E-mail: marmello@usp.br.

First published on April 25th, 2017 (English version).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4036687.svg)](https://doi.org/10.5281/zenodo.4036687)

Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again".

Disclaimer: You may use this software freely for any purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this software helps you produce any academic work (paper, book, chapter, dissertation, report, talk, lecture or similar), please acknowledge the authors and cite the source, using this repo's DOI or URL.


## List of folders and files

1. data (folder)

    a. net1.txt -> interactions of pollination between oil-collecting bees and oil flowers (Bezerra et al. 2009).

    b. net2.txt -> interactions of pollination between bats, hawkmoths, and plants (Queioz et al. 2020).

    c. net3an.txt -> interactions of seed destruction (antagonism) between maruspials, rodents, and plants (Genrich et al. 2017).
    
    d. net3mu.txt -> interactions of seed dispersal (mutualism) between maruspials, rodents, and plants (Genrich et al. 2017).
    
    
2. figures (folder)

    a. net1_bipartite_graph.png -> monolayer network drawn as a bipartite graph.

    b. net1_bipartite_matrix.png -> monolayer network drawn as an incidence matrix.
    
    c. net1_igraph_graph_multi.png -> multilayer network drawn as an energy-minimization graph, with node colors by taxonomic group.

    d. net1_igraph_graph_taxon.png -> monolayer network drawn as an energy-minimization graph, with node colors by taxonomic group.

    e. net1_igraph_graph.png -> monolayer network drawn as an energy-minimization graph.
    
    
3. network-drawing.pdf -> tutorial in PDF format.

4. network-drawing.R -> script with commented codes.

5. network-drawing.Rmd -> tutorial in RMD format.


## Data sources

1. net1 -> Bezerra, E. L. S., Machado, I. C., & Mello, M. A. R. (2009). [Pollination networks of oil-flowers: a tiny world within the smallest of all worlds](https://doi.org/10.1111/j.1365-2656.2009.01567.x). Journal of Animal Ecology, 78(5), 1096–1101. 

2. net2 -> Queiroz, J. A., Diniz, U. M., Vázquez, D. P., Quirino, Z. M., Santos, F. A. R., Mello, M. A. R., & Machado, I. C. (2020). Bats and hawkmoths form mixed modules with flowering plants in a nocturnal interaction network. Biotropica, *accepted*. See also this [repo](https://github.com/marmello77/queiroz-et-al-2020).

3. net3an and net3mu -> Genrich, C. M., Mello, M. A. R., Silveira, F. A. O., Bronstein, J. L., & Paglia, A. P. (2017). [Duality of interaction outcomes in a plant-frugivore multilayer network](https://doi.org/10.1111/oik.03825). Oikos, 126(3), 361–368. doi: 10.1111/oik.03825


## Instructions

If you are familiar with R, just run the scripts provided in "network-drawing.R".

If you are not so familiar with R, follow the steps described in the RMD tutorial.

If you just to take a look at the networks and their drawing solutions, read the PDF tutorial.


## Feedback

If you have any questions, corrections, or suggestions, please feel free to open and [issue](https://github.com/marmello77/network-drawing/issues) or make a [pull request](https://github.com/marmello77/network-drawing/pulls).


## Acknowledgments

We thank our labmates and our sponsors, especially the Alexander von Humboldt-Stiftung, CNPq, CAPES, and FAPESP, who gave us grants, fellowships, and scholarships. Pedro Jordano and Carsten Dormann helped us take our first steps in analyzing and drawing networks in R. Special thanks go to Katherine Ognyanova, who gave us invaluable tips on advanced network drawing in R. We strongly recommend her [online tutorial](http://kateto.net/network-visualization). Last, but not least, we thank the [Stack Overflow Community](https://stackoverflow.com), where we solve most of our coding dilemmas. 


## Suggested readings

* Ognyanova K. 2019. [Static and dynamic network visualization with R](http://kateto.net/network-visualization).

* Pocock, M. J. O., D. M. Evans, C. Fontaine, M. Harvey, R. Julliard, Ó. McLaughlin, J. Silvertown, A. Tamaddoni-Nezhad, P. C. L. White, and D. A. Bohan. (2016) [The Visualisation of Ecological Networks, and Their Use as a Tool for Engagement, Advocacy and Management](http://linkinghub.elsevier.com/retrieve/pii/S0065250415000355). In G. Woodward and D. A. Bohan (Eds.) Advances in Ecological Research. pp. 41–85, Academic Press, Cambridge.

* Marai, G. E., B. Pinaud, K. Bühler, A. Lex, and J. H. Morris. (2019) [Ten simple rules to create biological network figures for communication](https://doi.org/10.1371/journal.pcbi.1007244). PLOS Comput. Biol. 15: e1007244.

* Mello MAR, Muylaert RL, Pinheiro RBP & Félix GMF. 2016. [Guia para análise de redes ecológicas](https://marcomellolab.wordpress.com/software/). Edição dos autores, Belo Horizonte. 112 p. ISBN-13: 978-85-921757-0-2.

* Barabasi, A.L. (2016) [Network Science](http://barabasi.com/networksciencebook/), 1st ed. Cambridge University Press, Cambridge.

* Bascompte, J. & Jordano, P. (2014) [Mutualistic Networks](https://amzn.to/2FLwhto), 1st ed. Princeton University Press, Princeton.

