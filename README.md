<!-- README based on template: https://github.com/othneildrew/Best-README-Template -->

<a id="readme-top"></a>

<!-- PROJECT SHIELDS -->
[![Unlicense License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT HEADER -->
<br />
<div align="center">

  <h1 align="center">MSc internship repository - Tomer Sardjoe</h1>

  <p align="center">
    <i>Guanidine as the sole nitrogen source enables growth in _Alphaproteobacteria_ </i>
    <br />
    <a href=""><strong>Read the manuscript » (pending publication)</strong></a>
    <br />
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li>
      <a href="#usage">Usage</a>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- PROJECT DESCRIPTION -->
## Project description

This repository contains the raw data, and scripts that were written during my MSc internship project at the <a href="https://experts.exeter.ac.uk/25396-adam-monier"> Monier group</a> at <a href="https://www.exeter.ac.uk/research/institutes/livingsystems/"> the Living Systems Institute</a>. 

The main objective of the internship was to identifying whether guanidine hydrolase (gdmH) allows for the utilisation of guanidine as the sole nitrogen source in _Alphaproteobacteria_. The key findings of the project are:
* _gdmH_+ _Alphaproteobacteria_ grow on guanidine as the sole nitrogen source.
* _gdmH_ is widely, but sporadically distributed in _Alphaproteobacteria_.
* _gdmH_ occurs across diverse bacterial and eukaryotic lineages.


For more information and a closer look at the findings, please read the <a href=""> manuscript (pending publication)</a>.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

To be able to run the scripts and recreate the analyses, follow the steps below.

### Prerequisites

First, clone the repository onto your machine (note that the scripts were written and tested in a Linux environment, I will assume that you will clone the repository to a Linux system as well).
* Cloning the repository
```sh
git clone https://github.com/Tomersardjoe/MSc-internship.git
```

### Installation

The first thing you will want to do is to retrieve all the dependencies. For convenience, the project uses <a href="https://docs.conda.io/projects/conda/en/latest/index.html"> Conda</a> environments for this.

Create the Conda environments from the environments .yml files. These files can be found in the conda_envs directory and the environments can be created like so:
```sh
conda env create -f ncbi_tools.yml
conda env create -f phylo.yml
```
**The phylo environment should be activated for the R scripts, the ncbi_tools environment should be activated for all other scripts**.

Now, the project environment is ready for the analyses. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE -->
## Usage

The repository is divided into two main parts; data and scripts.

The data directory contains data that was generated from experiments in the lab or analyses based on those experiments (data that cannot be found anywhere else), and is therefore uploaded to the repository. Some of the data is used as input for analyses scripts.

The scripts directory contains support scripts, used to download or reformat data for subsequent analyses. The analyses directory contains the scripts that were used to generate the main figures in the manuscript.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GPL-3.0 license. See <a href="https://github.com/Tomersardjoe/MSc-internship/blob/main/LICENSE"> LICENSE</a> for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Tomer Sardjoe - tomer.sardjoe@wur.nl

Project Link: [https://github.com/Tomersardjoe/MSc-internship](https://github.com/Tomersardjoe/MSc-internship)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
The authors thank Marnix Medema for reviewing the manuscript and his helpful feedback. Additionally, we would like to express our gratitude to Victoria Jackson for transporting the Roscoff Culture Collection strains from France to the LSI.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/badge/License-GPL_3.0-green
[license-url]: https://github.com/Tomersardjoe/MSc-internship/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/LinkedIn-blue
[linkedin-url]: https://www.linkedin.com/in/tomersardjoe/
