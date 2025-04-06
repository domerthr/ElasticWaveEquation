# Elastic Wave Equation

This project provides an implementation of the **elastic wave equation** using space-time Galerkin finite elements. It supports the following discretization schemes:

- `cG(s)/cG(r)`
- `cG(s)/cG(r)-cG(r)`
- `cG(s)/dG(r)`

The goal is to provide these three methods and include support for initial and final time conditions for the velocity.

---

## How to Install

### 1. Install C++ Dependencies

To run the simulation, you need:

- **deal.II**  
  Installation instructions: [https://dealii.org/current/readme.html](https://dealii.org/current/readme.html)

- **ideal.II**  
  GitHub repository: [https://github.com/jpthiele/idealii](https://github.com/jpthiele/idealii)

### 2. Set Up Python Environment for Visualization

Python is used to visualize simulation results. You can install all necessary packages using the provided `environment.yml` file.

#### Create the Conda environment

```bash
conda env create -f environment.yml
```
---
## How to Use

### 1. Activate the Conda environment

```bash
conda activate elastic_wave
```
### 2. Clone and build the repository

```bash
git clone https://github.com/yourusername/elastic-wave-equation.git
```
