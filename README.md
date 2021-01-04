# Forecasting U.S 2020 Senate Election


<!-- TABLE OF CONTENTS -->
## Table of Contents

- [Forecasting U.S 2020 Senate Election](#forecasting-us-2020-senate-election)
  - [Table of Contents](#table-of-contents)
  - [About The Project](#about-the-project)
    - [Built With](#built-with)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Usage](#usage)
  - [License](#license)
  - [Contact](#contact)
  - [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project

This project builds a Dirichlet regression model for forecasting 2020 US Senate Election. We collect historical polling data and first build a Gaussian process model that captures the underlying time-evoling voter preference. The Dirichlet model then integrates the posteriors of voter preferences on election day and other election fundementals such as state-level partisan voting index, past experience of candidates and yearly swings. When forecasting 2020 races, our model naturally produces a posterior distribution of the future vote shares, rather than point estimates.


### Built With
* [Matlab](https://www.mathworks.com/products/matlab.html)
* [gpml3.6](http://gaussianprocess.org/gpml/code/matlab/release/oldcode.html)
  

<!-- GETTING STARTED -->
## Getting Started

To run this project, users need to have matlab and R installed locally. GPML matlab toolbox is already in the directory.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.
* gpml in matlab
```matlab
addpath("gpml-matlab-v3.6-2015-07-07");
startup;
```

### Installation

Clone the repo
```sh
git clone https://github.com/yunxiuqiu1115/CSE543T_FinalProject.git
```

<!-- USAGE EXAMPLES -->
## Usage

To forecast 2020 races at predefined horizons, run code
```sh
matlab -nodisplay -nodesktop -r "mainfunc(0); exit"
```
or run command below in matlab.
```sh
mainfunc(0)
```

The prediction plots of each candidates can be found in `plots` repository and the posttrain results are in `results` repository.


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Yehu Chen - chenyehu@wustl.edu

Yunxiu Qiu - yunxiuqiu@wustl.edu

Yanpeng Yuan - yanpeng@wustl.edu

Project Link: [https://github.com/yunxiuqiu1115/CSE543T_FinalProject/tree/dev](https://github.com/yunxiuqiu1115/CSE543T_FinalProject/tree/dev)


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [gpml](http://www.gaussianprocess.org/gpml/code/matlab/doc/)
* [GitHub Pages](https://pages.github.com)

