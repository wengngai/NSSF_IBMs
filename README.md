# NSSF_IBMs
This project aims to construct spatially explicit individual based models of tree species in the Nee Soon Freshwater Swamp Forest (NSSF) catchment in Singapore, and to use these models to predict tree community compositions in a heterogeneous freshwater swamp forest landscape.
Models are constructed from vital rates estimated in an earlier study on habitat mismatch (Citation here).

## Species
1. Prunus polystachya (Rosaceae) a non-swamp specialist
2. Strombosia ceylanica (Olacaceae) a generalist
3. Pometia pinnata (Sapindaceae) a swamp specialist

(Gironniera nervosa (Cannabaceae) was originally considered but because it is dioecious, it was dropped.)

## Tree life stages
Each tree passes through two life stages: a seedling stage, where the unit of size is height (cm; log-transformed), and an adult stage, where the unit of size is DBH (cm; log-transformed).

## Vital rate and other models
These are predictors from the top models for each vital rate, as reported in (earlier study citation here)
1. Seedling (DBH < 1 cm) growth: conspecific and heterospecific adult effects, and their interactions with seedling height
2. Small tree (seedlings + saplings < 2.7 cm DBH) survival*: intraspecific competition and habitat mismatch
3. Adult growth**: intra- and interspecific competition
4. Large tree (DBH > 2.7 cm) survival*: intra- and interspecific competition
5. Fruiting probability: habitat mismatch

*Small and large tree survival models are based off models from Needham et al. (Proc. B. 2018) and Johnson et al. (Nat. Ecol. Evol. 2018)
**Adult growth models are based off the model of Kohyama et al. (J. Ecol. 2020)
These models are presented in the RShiny app file

Several other models also provide the foundation for constructing these IBMs: 
1. Seedling-sapling allometry model. This model is used for transitions between seedling and adult phases, and maps seedling log-height to adult log-DBH
2. Dispersal kernels. This model was constructed from seedling data from belt transects from isolated parent trees, and allows the estimation of locations of recruitment events around reproducing trees. The process of estimating recruitment rates and locations are described in the file (insert name here)
