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
    <a href="https://github.com/shuyuzheng/TidyComb"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://drugcomb.fimm.fi">DrugComb Portal</a>
    ·
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
install_github("shuyuzheng/TidyComb")
```

## Usage


## Roadmap

## Contributing

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
