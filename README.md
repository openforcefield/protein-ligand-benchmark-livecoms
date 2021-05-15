# Best practices for constructing, preparing, and evaluating protein-ligand binding affinity benchmarks

## Abstract

Free energy calculations are rapidly becoming indispensable in structure-enabled drug discovery programs. As new methods, force fields, and implementations are developed, assessing their expected accuracy on real-world systems (benchmarking) becomes critical to provide users with an assessment of the accuracy expected when these methods are applied within their domain of applicability, and developers with a way to assess the expected impact of new methodologies. These assessments require construction of a benchmark—a set of well-prepared, high quality systems with corresponding experimental measurements designed to ensure the resulting calculations provide a realistic assessment of expected performance when these methods are deployed within their domains of applicability. To date, the community has not yet adopted a common standardized benchmark, and existing benchmark reports suffer from a myriad of issues, including poor data quality, limited statistical power, and statistically deficient analyses, all of which can conspire to produce benchmarks that are poorly predictive of real-world performance. Here, we address these issues by presenting guidelines for 

1. curating experimental data to develop meaningful benchmark sets, 
2. preparing benchmark inputs according to best practices to facilitate widespread adoption, and 
3. analysis of the resulting predictions to enable statistically meaningful comparisons among methods and force fields. 
 
We highlight challenges and open questions that remain to be solved in these areas, as well as recommendations for the collection of new datasets that might optimally serve to measure progress as methods become systematically more reliable. Finally, we provide a curated, versioned, open, standardized benchmark set adherent to these standards [protein-ligand-benchmark](https://github.com/openforcefield/protein-ligand-benchmark)  and an open source toolkit for implementing standardized best practices assessments [arsenic](https://github.com/openforcefield/Arsenic) for the community to use as a standardized assessment tool. While our main focus is free energy methods based on molecular simulations, these guidelines should prove useful for assessment of the rapidly growing field of machine learning methods for affinity prediction as well.

## List of Authors

DH: David Hahn    
CIB: Christopher I. Bayly  
HBM: Hannah Bruce Macdonald
JDC: John D. Chodera
ASJSM: Antonia S. J. S. Mey
DLM:  David L. Mobley
LPB: Laura Perez Benito
CS Christina E.M. Schindler
GT: Gary Tresadern
GW: Gregory L. Warren

## List of Contributors
(non-author list of people who contributed to document)

## Paper Writing as Code Development
<!-- This discussion is so that people know how to contribute to your document. -->
This paper is being developed as a living document, open to changes from the community.
You can read more about the concept of writing a paper in the same way one would write software code in the essay ["Paper writing as code development"](https://livecomsjournal.github.io/paper_code.html).
If you have comments or suggestions, we welcome them! Please [submit them as issues](https://guides.github.com/features/issues/) to this GitHub repository so they can be recorded and given credit for the contribution.
Specific changes can be proposed [via pull requests](https://help.github.com/articles/about-pull-requests/).


## Online Resources
Original brainstorming document: (https://docs.google.com/document/d/1lCGcol6jYLQmcfqrUv9h_FsWygTZzqYxqgjOLCyMoL4/edit)

## Design Notes:
- Figure labels are in captital bold face
- Standard font for all figures should be Helvetica. Font size for lables should be at least the size of the document text. 

## List of Released Versions
still working on first release
