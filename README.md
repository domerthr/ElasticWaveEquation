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

  For the installation of `ideal.II`, make sure to include the following adjustments:

  - Add **1D case** support in all `*.inst` files.
  - Modify `spacetime_dof_handler.hh` and `spacetime_dof_handler.cc` by adding the `void reinit()` function for handling degrees of freedom:

    ```cpp
    template <int dim>
    void
    DoFHandler<dim>::reinit()
    {
      if (_tria != nullptr)
      {
        this->_dof_handlers.clear();
        slab::TriaIterator<dim> tria_it = this->_tria->begin();
        slab::TriaIterator<dim> tria_end = this->_tria->end();
        for (; tria_it != tria_end; ++tria_it)
        {
          this->_dof_handlers.push_back(
            idealii::slab::DoFHandler<dim>(*tria_it));
        }
      }
      else
      {
        Assert(false, dealii::ExcInternalError());
      }
    }
    ```
Make sure you have all features installed (parallel linear algebra and other advanced features)

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
git clone https://github.com/domerthr/ElasticWaveEquation.git
cmake -S. -Bbuild -DIDEAL_II_DIR=<path_to_install_idealii_in>
```
### 3. Run a simulation
```bash
cd build
make run
```
