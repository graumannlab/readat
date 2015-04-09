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
