[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/0.1.0/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://semaphoreci.com/api/v1/richierocks/readat/branches/master/badge.svg)](https://semaphoreci.com/richierocks/readat)
[![Build status](https://ci.appveyor.com/api/projects/status/7lbjxa2snrhbgjcr?svg=true)](https://ci.appveyor.com/project/richierocks/readat)

[//]: # ([Bioconductor build status](http://bioconductor.org/shields/build/release/bioc/Biobase.svg))
[//]: # ([Time on Bioconductor](http://bioconductor.org/shields/years-in-bioc/BiocGenerics.svg))

# readat

Tools for importing and working with [SomaLogic](http://www.somalogic.com/Homepage.aspx) [ADAT](https://bitbucket.org/graumannlabtools/adat-spec) files.

Read the paper in BMC Bioinformatics: [readat: An R package for reading and working with SomaLogic ADAT files](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1007-8)

### Installation

To install the package, you first need the 
[*devtools*](https://github.com/hadley/devtools) package.

```{r}
install.packages("devtools")
```

Then you can install the *readat* package using

```{r}
devtools::install_bitbucket("graumannlabtools/readat")
```

# Functionality

`readAdat` let's you import data from the SomaLogic ADAT file format.  The
result is stored in an object of class `WideSomaData`.  This is a data.table 
with one sample per row, and includes both the intensities and sample metadata.
The object also has an attribute named `sequenceData` that contains a data.table
of sequence metadata.  There are also attributes for experiment metadata and a 
checksum.

`getSequenceData`, `getMetadata`, and `getChecksum` provide shortcuts to access 
these attributes.

There is a `melt` method for `WideSomaData` that converts the object to 
`LongSomaData` format, which is a data.table with one intensity per row.

`soma2eset` converts an object of class `WideSomaData` to a 
`Biobase::ExpressionSet`.

You can retrieve additional annotations by SomaLogic Sequence ID uisng 
`getEnsemblIds`, `getUniProtKeywords`, `getChromosomalPositions`, `getPfam`, 
`getKeggDefinitions`, `getKeggModules`, `getKeggPathways`, 
`getGoMolecularFunctions`, `getGoBiologicalProcesses`, `getGoCellularComponents`.

`getSequencesWithLargestBetweenGroupVariation` find the sequences with the 
largest amount of variation between specified sample groups.
