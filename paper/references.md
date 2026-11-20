# References for the manuscript

Each entry was looked up on 2026-11-18 (PubMed unless marked otherwise). Check the final list against the journal's style before submission. "Use" says where the reference belongs in the paper.

## Software the tool builds on

| Reference | PMID | DOI | Use |
|---|---|---|---|
| Schindelin J et al. Fiji: an open-source platform for biological-image analysis. Nat Methods 2012;9:676-682. | 22743772 | 10.1038/nmeth.2019 | Platform |
| Schneider CA, Rasband WS, Eliceiri KW. NIH Image to ImageJ: 25 years of image analysis. Nat Methods 2012;9:671-675. | 22930834 | 10.1038/nmeth.2089 | Platform |
| Linkert M et al. Metadata matters: access to image data in the real world. J Cell Biol 2010;189:777-782. | 20513764 | 10.1083/jcb.201004104 | Bio-Formats |
| Phansalkar N, More S, Sabale A, Joshi M. Adaptive local thresholding for detection of nuclei in diversity stained cytology images. ICCSP 2011. | (not in PubMed; web-verified) | 10.1109/ICCSP.2011.5739305 | Local threshold |
| Zuiderveld K. Contrast limited adaptive histogram equalization. In: Graphics Gems IV, 1994, pp. 474-485. | (book chapter; web-verified, no DOI) | | CLAHE |
| Legland D, Arganda-Carreras I, Andrey P. MorphoLibJ. Bioinformatics 2016;32:3532-3534. | 27412086 | 10.1093/bioinformatics/btw413 | Distance maps (related) |

## Comparator and related tools

| Reference | PMID | DOI | Use |
|---|---|---|---|
| Schmidt U, Weigert M, Broaddus C, Myers G. Cell detection with star-convex polygons. MICCAI 2018, LNCS 11071:265-273. | (web-verified) | 10.1007/978-3-030-00934-2_30 | Comparator (StarDist) |
| Stringer C et al. Cellpose: a generalist algorithm for cellular segmentation. Nat Methods 2021;18:100-106. | 33318659 | 10.1038/s41592-020-01018-x | Comparator |
| Pachitariu M, Stringer C. Cellpose 2.0. Nat Methods 2022;19:1634-1641. | 36344832 | 10.1038/s41592-022-01663-4 | Comparator |
| Stirling DR et al. CellProfiler 4. BMC Bioinformatics 2021;22:433. | 34507520 | 10.1186/s12859-021-04344-9 | State of the field |
| McQuin C et al. CellProfiler 3.0. PLoS Biol 2018;16:e2005970. | 29969450 | 10.1371/journal.pbio.2005970 | State of the field |
| Bankhead P et al. QuPath. Sci Rep 2017;7:16878. | 29203879 | 10.1038/s41598-017-17204-5 | State of the field |
| Berg S et al. ilastik. Nat Methods 2019;16:1226-1232. | 31570887 | 10.1038/s41592-019-0582-9 | State of the field |
| Beretta CA et al. Quanty-cFOS. Cells 2023;12:704. | 36899840 | 10.3390/cells12050704 | Closest neighbour, validation precedent |
| Unger et al. SAM. Sci Rep 2025. | 40721622 | 10.1038/s41598-025-12144-x | Closest neighbour |
| Waisman. Multi-Channel Nuclear Analysis. MicroPubl Biol 2025. | 40547381 | 10.17912/micropub.biology.001590 | Closest neighbour |
| Sanchez-Arias JC et al. Front Neuroanat 2021. | 34949993 | 10.3389/fnana.2021.722443 | Best practice, cell density in sections |
| Courtney JM et al. eNeuro 2021. | 34642225 | 10.1523/ENEURO.0177-21.2021 | Threshold sweep against manual counts |
| Schroder et al. STAR Protoc 2024. | 38427568 | 10.1016/j.xpro.2024.102904 | Organoid-section protocol precedent |
| Gritti N et al. MOrgAna. Development 2021. | 34494114 | 10.1242/dev.199611 | Organoid tools |
| Mukashyaka P et al. Cellos. Nat Commun 2023. | 38114489 | 10.1038/s41467-023-44162-6 | Organoid tools |

## Prior art of each module

| Module | Reference | PMID | DOI |
|---|---|---|---|
| Depth profile | Albanese A et al. Multiscale 3D phenotyping of human cerebral organoids (SCOUT). Sci Rep 2020;10:21487. | 33293587 | 10.1038/s41598-020-78130-7 |
| Depth profile | Spoerri L et al. Front Digit Health 2021. | 34713141 | 10.3389/fdgth.2021.668390 |
| Depth profile | Ciarpella F et al. iScience 2021. | 34901791 | 10.1016/j.isci.2021.103438 |
| Threshold spectrum | Johnson SJ, Walker FR. Strategies to improve quantitative assessment of immunohistochemical and immunofluorescent labelling. Sci Rep 2015;5:10607. | 26066532 | 10.1038/srep10607 |
| Threshold spectrum | Choudhury KR et al. J Histochem Cytochem 2010. | 19687472 | 10.1369/jhc.2009.953554 |
| Auto-tune | Held C et al. Approaches to automatic parameter fitting in a microscopy image segmentation pipeline. J Pathol Inform 2013;4:S5. | 23766941 | 10.4103/2153-3539.109831 |
| Auto-tune | Teodoro G et al. Algorithm sensitivity analysis and parameter tuning for tissue image segmentation pipelines. Bioinformatics 2017;33:1064-1072. | 28062445 | 10.1093/bioinformatics/btw749 |
| Auto-tune | Taveira LFR et al. J Digit Imaging 2019;32:521-533. | 30402669 | 10.1007/s10278-018-0138-z |
| Corrected counts | Ritch MD et al. AxoNet. Sci Rep 2020;10:8034. | 32415269 | 10.1038/s41598-020-64898-1 |
| Corrected counts | Reynaud J et al. Invest Ophthalmol Vis Sci 2012;53:2951-2959. | 22467571 | 10.1167/iovs.11-9274 |
| Corrected counts | Angelopoulos AN et al. Prediction-powered inference. Science 2023. | 37943906 | 10.1126/science.adi6000 |
| Robustness | Steegen S et al. Increasing transparency through a multiverse analysis. Perspect Psychol Sci 2016. | 27694465 | 10.1177/1745691616658637 |
| Robustness | Simonsohn U, Simmons JP, Nelson LD. Specification curve analysis. Nat Hum Behav 2020. | 32719546 | 10.1038/s41562-020-0912-z |
| Robustness | Botvinik-Nezer R et al. Variability in the analysis of a single neuroimaging dataset by many teams. Nature 2020. | 32483374 | 10.1038/s41586-020-2314-9 |
| Robustness | Gomes et al. Comput Biol Med 2019 (segmentation parameters change survival conclusions). | 31054503 | 10.1016/j.compbiomed.2019.03.006 |

## Validation and reporting standards

| Reference | PMID | DOI |
|---|---|---|
| Maier-Hein L et al. Metrics reloaded. Nat Methods 2024. | 38347141 | 10.1038/s41592-023-02151-z |
| Reinke A et al. Understanding metric-related pitfalls. Nat Methods 2024. | 38347140 | 10.1038/s41592-023-02150-0 |
| Caicedo JC et al. Evaluation of deep learning strategies for nucleus segmentation. Cytometry A 2019. | 31313519 | 10.1002/cyto.a.23863 |
| Bland JM, Altman DG. Statistical methods for assessing agreement. Lancet 1986. | 2868172 | (none returned) |
| Koo TK, Li MY. ICC guideline. J Chiropr Med 2016. | 27330520 | 10.1016/j.jcm.2016.02.012 |
| Schmied C et al. Community-developed checklists for publishing images and image analyses. Nat Methods 2024. | 37710020 | 10.1038/s41592-023-01987-9 |
| Abercrombie M. Estimation of nuclear population from microtome sections. Anat Rec 1946. | 21015608 | 10.1002/ar.1090940210 |

## Not verified

- Requirements of SoftwareX, BMC Bioinformatics and STAR Protocols: their author pages could not be opened. Read them directly before choosing.
- The statement that the area under the full coverage curve equals the normalized mean intensity is a derivation (the integral of a survival function), not a citation. It is checked by a test in `tests/test_spectrum.py`.
