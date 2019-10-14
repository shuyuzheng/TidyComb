[![Contributors][contributors-shield]][contributors-url]
[![MPL2.0 License][license-shield]][license-url]
![Twitter Follow](https://img.shields.io/twitter/follow/DrugComb.svg?style=social)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <h3 align="center">TidyComb</h3>

  <p align="center">
    Tidy drug combination high-throughput screening data for DrugComb portal
    <br />
    <a href="https://drugcomb.fimm.fi"><strong>Explore DrugComb Portal »</strong></a>
    <br />
    <a href="https://github.com/shuyuzheng/TidyComb/issues">Report Bug</a>
    ·
    <a href="https://github.com/shuyuzheng/TidyComb/issues">Request Feature</a>
  </p>
</p>


<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Package](#about-the-package)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Reference](#reference)
* [License](#license)
* [Contact](#contact)

## About the Package 

The R package **TidyComb** wrapped the functions that could be used to tidy drug combination data from various sources into the format that is required by uploading to [DrugComb](https://drugcomb.fimm.fi) database. It includes functions for retrieving drug information from chemistry databases: PubChem and Chembl, retrieving cell line information from Cellosaurus, calculating drug synergy scores, etc.

### Built with

* [R](https://www.r-project.org/) 3.6
* [RStudio](https://www.rstudio.com/) 1.2

## Getting Start

### Prerequistes

**TidyComb** is a package written with [R](https://www.r-project.org/) programming language. Please make sure that you have installed [R project](https://www.r-project.org/) ( >= V3.6) on your local machine. You can download it from [CRAN](https://cran.r-project.org/mirrors.html).

### Installation

You can install **TidyComb** from this [GitHub repository](https://github.com/shuyuzheng/TidyComb):

1. Install and load the [devtools](https://github.com/hadley/devtools) package. You can do this from [CRAN](https://cran.r-project.org/). Invoke R and then type

```
install.packages("devtools")
library(devtools)
```

2. To install **TidyComb** from [GitHub](https://github.com/), you'd type:

```
install_github("DrugComb/TidyComb")
```

## Usage

The functions in **TidyComb** can be devided into 2 categories:

1. Functions for collecting information from external databases:

* [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed)
* [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
* [ChEMBL](https://www.ebi.ac.uk/chembl/)
* [UniChem](https://www.ebi.ac.uk/unichem/)
* [DrugBank](https://www.drugbank.ca/)
* [Cellosaurus](https://web.expasy.org/cgi-bin/cellosaurus/search)

The APIs used in those functions:

* [PubMed Central APIs](https://www.ncbi.nlm.nih.gov/pmc/tools/developers/)
* [PubChem PUG REST](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
* [ChEMBL web services API](https://www.ebi.ac.uk/chembl/api/data/docs)
* [UniChem RESTful Web Service API](https://www.ebi.ac.uk/unichem/info/webservices)

2. Functions for analysizing drug synergy and sensitivity. Synergy scores are calculated using four different reference models:

* Bliss model [Bliss, 1939][2] assumes a stochastic process in which two drugs elicit their effects independently, and the expected combination effect can be calculated based on the probability of independent events 
* Highest Single Agent (HSA) [Berenbaum, 1989][3] states that the expected combination effect equals to the higher effect of individual drugs
* Loewe additivity model [Loewe, 1953][4] defines the expected effect yLOEW E as if a drug was combined with itself. Unlike the HSA and the Bliss independence models giving a point estimate using different assumptions, the Loewe additivity model considers the dose-response curves of individual drugs.
* Zero Interaction Potency (ZIP) [Yadav et al., 2015][5] calculates the expected effect of two drugs under the assumption that they do not potentiate each other.

Sensitivity is calculated by CSS [Malyutina et al., 2019][6] - drug combination sensitivity score is derived using relative IC50 values of compounds and the area under their dose-response curves. 

## Reference

[1]: Zagidullin B., Aldahdooh J., Zheng S., Wang W., Wang Y., Saad J., Malyutina A., Jafari M., Tanoli Z., Pessia A., Tang J. (2019). DrugComb: an integrative cancer drug combination data portal, Nucleic Acids Research, 47(W1):W43-W51.
[2]: Bliss, C. I. (1939). The toxicity of poisons applied jointly1. Annals of Applied Biology, 26(3):585–615.
[3]: Berenbaum, M. C. (1989). What is synergy? Pharmacol. Rev., 41(2):93–141.
[4]: Loewe, S. (1953). The problem of synergism and antagonism of combined drugs. Arzneimit- telforschung, 3(6):285–290.
[5]: Yadav, B., Wennerberg, K., Aittokallio, T., and Tang, J. (2015). Searching for Drug Synergy in Complex Dose-Response Landscapes Using an Interaction Potency Model. Comput Struct Biotechnol J, 13:504– 513.
[6]: Malyutina A., Majumder MM., Wang W., Pessia A., Heckman CA., Tang J. (2019). Drug combination sensitivity scoring facilitates the discovery of synergistic and efficacious drug combinations in cancer. PLoS Comput Biol, 15(5):e1006752. 

## License

Distributed under the Mozilla Public License 2.0

## Contact

Shuyu Zheng - shuyu.zheng@helsinki.fi

Project Link: https://github.com/shuyuzheng/TidyComb

## Acknowledgements
* [Img Shields](https://shields.io)
* [Choose an Open Source License](https://choosealicense.com)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/badge/contributors-1-orange.svg?style=flat-square
[contributors-url]: https://github.com/shuyuzheng/TidyComb/graphs/contributors
[license-shield]: https://img.shields.io/badge/license-MPL--2.0-blue.svg
[license-url]: https://choosealicense.com/licenses/mpl-2.0
