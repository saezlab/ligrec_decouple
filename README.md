# Systematic Comparison of Cell-Cell Communication Tools and Resources
  
![landingpage](ligrec_pipe.png)

## Objectives
  
The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked
an immense interest in understanding intercellular crosstalk. Multiple
tools and resources that aid the investigation of cell-cell communication (CCC)
were published recently.
However, these methods and resources are usually in a fixed combination of a
tool and its corresponding resource, but in principle any resource could be
combined with any method. Yet, it is largely unclear the
difference that the choice of resource and tool can have on the predicted
CCC events. Thus, we attempt to shed some light on this topic via the
systematic comparison of how different combinations might influence CCC
inference.
  
  

## LIANA framework
  
To this end we built [LIANA](https://github.com/saezlab/liana), a framework to decouple the tools from their corresponding resources.
  
  
## Tools

The tools implemented in this repository are:

- CellPhoneDB algorithm (via [Squidpy](https://squidpy.readthedocs.io/en/latest/))
- CellChat
- NATMI
- Connectome
- SingleCellSignalR (SCA)
- iTALK
  
  
## Resources

The following CCC resources are accessible via this pipeline:

- CellChatDB
- CellPhoneDB
- Ramilowski2015
- Baccin2019
- LRdb
- Kiroauc2010
- ICELLNET
- iTALK
- EMBRACE
- HPMR
- Guide2Pharma
- connectomeDB2020
- talklr
- CellTalkDB
- OmniPath
  
  
### OmniPath
  
All the resources above are retrieved from OmniPath (https://omnipathdb.org/).
OmniPath itself is also a composite resource combining all the ones listed
above. However the cell-cell interactions in OmniPath are more than simply
the union of the ligand-receptor resources. OmniPath uses several further
databases to collect information about the roles of proteins in intercellular
communication and again other databases to find connections between them. At
the same time, OmniPath blacklists certain wrong annotations, removing some
of the contents of the original resources. However the data of individual
resources retrieved from the OmniPath web service is not affected by this,
each resource supposed to be identical to its original form, apart from minor
processing imperfections.
  
  
### Reshuffled and Default Resources
  
Moreover, a Reshuffled resource can be generated via reshuffling any of the
abovementioned using the `BiRewire` package, and each tool can be run with
its 'Default' resource, the dataset used in its original publication.
  
  

   
   
```{r}
library(devtools)  
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK 
install_github('saezlab/OmnipathR') # install latest version of OmnipathR
devtools::install_github("sqjin/CellChat")  
devtools::install_github('msraredon/Connectome', ref = 'master')   
devtools::install_github('Coolgenome/iTALK', ref = 'biorxiv')   
# A modified version of SingleCellSignalR (SCA) that enables external resources
devtools::install_github(repo = "CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")
# Install this repos
devtools::install_github('saezlab/ligrec_decoupler')
```
  
  
#### Clone NATMI, as it is run via the command line 
```{sh}
git clone https://github.com/asrhou/NATMI.git
```
  
  
## References
Baccin, C., Al-Sabah, J., Velten, L., Helbling, P.M., Grünschläger, F., Hernández-Malmierca, P., Nombela-Arrieta, C., Steinmetz, L.M., Trumpp, A., and Haas, S. (2020). Combined single-cell and spatial transcriptomics reveal the molecular, cellular and spatial bone marrow niche organization. Nat. Cell Biol. 22, 38–48.

Ben-Shlomo, I., Yu Hsu, S., Rauch, R., Kowalski, H.W., and Hsueh, A.J.W. (2003). Signaling receptome: a genomic and evolutionary perspective of plasma membrane receptors involved in signal transduction. Sci. STKE 2003, RE9.

Cabello-Aguilar, S., Alame, M., Kon-Sun-Tack, F., Fau, C., Lacroix, M., and Colinge, J. (2020). SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics. Nucleic Acids Res. 48, e55.

Efremova, M., Vento-Tormo, M., Teichmann, S.A., and Vento-Tormo, R. (2020). CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat. Protoc. 15, 1484–1506.

Harding, S.D., Sharman, J.L., Faccenda, E., Southan, C., Pawson, A.J., Ireland, S., Gray, A.J.G., Bruce, L., Alexander, S.P.H., Anderton, S., et al. (2018). The IUPHAR/BPS Guide to PHARMACOLOGY in 2018: updates and expansion to encompass the new guide to IMMUNOPHARMACOLOGY. Nucleic Acids Res. 46, D1091–D1106.

Hou, R., Denisenko, E., Ong, H.T., Ramilowski, J.A., and Forrest, A.R.R. (2020). Predicting cell-to-cell communication networks using NATMI. Nat. Commun. 11, 5011.

Jin, S., Guerrero-Juarez, C.F., Zhang, L., Chang, I., Ramos, R., Kuan, C.-H., Myung, P., Plikus, M.V., and Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 12, 1088.


Noël, F., Massenet-Regad, L., Carmi-Levy, I., Cappuccio, A., Grandclaudon, M., Trichot, C., Kieffer, Y., Mechta-Grigoriou, F., and Soumelis, V. (2021). Dissection of intercellular communication using the transcriptome-based framework ICELLNET. Nat. Commun. 12, 1089.

Palla, G., Spitzer, H., Klein, M., Fischer, D.S., Schaar, A.C., Kuemmerle, L.B., Rybakov, S., Ibarra, I.L., Holmberg, O., Virshup, I., et al. (2021). Squidpy: a scalable framework for spatial single cell analysis. BioRxiv.

Ramilowski, J.A., Goldberg, T., Harshbarger, J., Kloppmann, E., Lizio, M., Satagopam, V.P., Itoh, M., Kawaji, H., Carninci, P., Rost, B., et al. (2015). A draft network of ligand-receptor-mediated multicellular signalling in human. Nat. Commun. 6, 7866.

Raredon, M.S.B., Yang, J., Garritano, J., Wang, M., Kushnir, D., Schupp, J.C., Adams, T.S., Greaney, A.M., Leiby, K.L., Kaminski, N., et al. (2021). Connectome: computation and visualization of cell-cell signaling topologies in single-cell systems data. BioRxiv.

Shao, X., Liao, J., Li, C., Lu, X., Cheng, J., and Fan, X. (2020). CellTalkDB: a manually curated database of ligand-receptor interactions in humans and mice. Brief. Bioinformatics.

Sheikh, B.N., Bondareva, O., Guhathakurta, S., Tsang, T.H., Sikora, K., Aizarani, N., Sagar, Holz, H., Grün, D., Hein, L., et al. (2019). Systematic Identification of Cell-Cell Communication Networks in the Developing Brain. IScience 21, 273–287.

Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., et al. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Mol. Syst. Biol. 17, e9923.

Vento-Tormo, R., Efremova, M., Botting, R.A., Turco, M.Y., Vento-Tormo, M., Meyer, K.B., Park, J.-E., Stephenson, E., Polański, K., Goncalves, A., et al. (2018). Single-cell reconstruction of the early maternal-fetal interface in humans. Nature 563, 347–353.

Wang, Y. (2020). talklr uncovers ligand-receptor mediated intercellular crosstalk. BioRxiv.

Wang, Y., Wang, R., Zhang, S., Song, S., Jiang, C., Han, G., Wang, M., Ajani, J., Futreal, A., and Wang, L. (2019). iTALK: an R Package to Characterize and Illustrate Intercellular Communication. BioRxiv.

