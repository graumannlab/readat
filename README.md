# somalogic

Tools for importing and working with SomaLogic ADAT files.

### Installation

To install the package, you first need the 
[*devtools*](https://github.com/hadley/devtools) package.

```{r}
install.packages("devtools")
```

Then you can install the *somalogic* package using

```{r}
library(devtools)
install_bitbucket(
  "graumannlab/somalogic",
  auth_user = "your bitbucket username", 
  password  = "your bitbucket password"  
)
```

# Functionality

`readSomaLogic` let's you import data from the SomaLogic ADAT file format.  The
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




