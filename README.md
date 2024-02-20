# Synthetic thermal turbulent inflow generation (Extended SEM)

## This version is v0.9-beta. v1.0 will be announced later. 

The extended Synthetic-Eddy Method (XSEM) offers an efficient and robust computational procedure for generating realistic thermal turbulent inflow for eddy-resolving simulations of inflow-outflow problems. Initially proposed by Jarrin in 2006, this method was extended by Oh et al. in 2019 to include thermal field fluctuations alongside turbulent ones. The XSEM combines coherent structures with prescribed mean fields, produced through an eddy-convecting process within a virtual box. This process ensures alignment of the velocity component and temperature fluctuations with the second-order statistics of turbulent thermal boundary layers.

## Authors

- Geunwoo Oh (gwoh@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University.
- Kyung Min Noh (blackcata109@gmail.com), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University.
- Hyunwook Park (deukgyun@naver.com), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University.
- Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University.
- Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University.

## Usage
### Downloading PaScaL_TDMA
The repository can be cloned as follows:

```
git clone https://github.com/MPMC-Lab/PaScaL_TDMA.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.


### Compile and Build

#### Prerequisites

- A Fortran compiler is required to compile XSEM.

#### Build Instructions

1. **Build XSEM**:
   ```
   make
   ```
2. **Build an example problem after building XSEM**:
   ```
   make example
   ```

### Running the Example

After building, an executable binary `a.out` is created in the `run` folder. Use the pre-defined input file `PARA_INPUT.inp` in the run folder to execute the example:
```
./a.out
```

### Post-processing

The results from executing `a.out` will be saved in the `run` folder.

## Folder Structure

- `Src`: Source files of XSEM.
- `Example`: Source files for an example problem related to urban-scale boundary layers.
- `Run`: Contains the executable binary file for the example problem after building.
- `Doc`: Documentation.

## Citation

Please cite the following papers if you refer to this project:

```bibtex
@article{oh2019extended,
  title={Extended synthetic eddy method to generate inflow data for turbulent thermal boundary layer},
  author={Oh, Geunwoo and Noh, Kyung Min and Park, Hyunwook and Choi, Jung-Il},
  journal={International Journal of Heat and Mass Transfer},
  volume={134},
  pages={1261--1267},
  year={2019},
  publisher={Elsevier}
}

@article{yang2023multi,
  title={Multi-GPU-based real-time large-eddy simulations for urban microclimate},
  author={Yang, Mingyu and Oh, Geunwoo and Xu, Tiantian and Kim, Jungwoo and Kang, Ji-Hoon and Choi, Jung-Il},
  journal={Building and Environment},
  volume={245},
  pages={110856},
  year={2023},
  publisher={Elsevier}
}
```

## Reference

For more information, please refer to the reference paper and the Multi-Physics Modeling and Computation (MPMC) Lab.
