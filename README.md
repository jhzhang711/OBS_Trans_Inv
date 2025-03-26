# OBS_Trans_Inv
Seismological Transdimensional Bayesian Inversion Tools for OBS in Marine Environments.

## 📌 

This is a transdimensional Bayesian inversion tool applicable to marine environments. 
The code largely follows the methodology described in Bodin et al. (2012). However, it additionally accounts for the influence of the overlying water layer on ocean bottom seismometers(OBS) by incorporating the water layer thickness as a prior in the inversion process. We have referred to the work of Akuhara et al. (2023) in this section.

## 📦 

```

├── codes-inv/         # codes
├── RFx.obs            # RF Observations
├── SWD_xxxxx.obs      # Dispersion  Observations
├── REF_in.mod         # Initial Model for Inversion
├── inputmodel.txt     # Model Used to Generate Observations
├── plotResult.sh      # Script to Plot Results 
├── .gitignore         # Git 忽略文件
├── README.md          # 说明文档
└── LICENSE            # 许可证
```

## 🚀 Installation 
Type make in the codes-inv directory. Please edit the Makefile in accordance with your environment (i.e., compiler type and libarary paths).

### 🔧 Requirements

```sh
# for example
module load mpich/3.1.4-gcc4.9.2  
module load fftw/3.3.8-mpi
module load lapack
```

### 🏃 How to run

```sh
# for example
mkdir -p ./OUT
mpirun ./codes-inv/JointINV > ${STA}.out
```

## 📜 Quick Guidance

You can modify the prior information of the water layer in Joint.f90, as well as specify the inversion type, number of data points, weights, iterations, and other parameters. 
Additionally, you can adjust the prior velocity and velocity ratio information in priorvalue.f90 and priorvpvs.f90.

## 📄 

GNU General Public License 


