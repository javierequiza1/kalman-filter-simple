# High-Performance Kalman Filter & Hyperparameter Optimization 🚀

![Fortran](https://img.shields.io/badge/Fortran-Modern-734f96?style=flat-square&logo=fortran)
![Gnuplot](https://img.shields.io/badge/Visualization-Gnuplot-red?style=flat-square)
![Status](https://img.shields.io/badge/Status-Optimization_Complete-success?style=flat-square)

A robust implementation of a **1D Kalman Filter** in Modern Fortran designed to track oscillating signals with high precision. 

This project goes beyond standard filtering by implementing a **Grid Search Optimizer** that generates a 3D error hypersurface. This allows for the precise identification of the optimal Process Noise Covariance matrices ($Q$) by minimizing the Root Mean Square Error (RMSE).

## 🌟 Key Features

* **⚡ High Performance:** Written in optimized Modern Fortran using vector operations.
* **🧠 Branchless Logic:** Utilizes Fortran's `merge` intrinsic for conditional logic, ensuring constant execution time and pipeline efficiency (avoiding `IF` branching in the main loop).
* **🎯 Hyperparameter Optimization:** Includes a module to explore the hyperparameter space ($\log_{10} Q_{pos}$ vs $\log_{10} Q_{vel}$), generating a heat map of the error landscape.
* **📊 Professional Visualization:** Scripts for **Gnuplot** to render real-time tracking dashboards and 3D optimization surfaces.

## 📈 Visualizations

### 1. The Optimization Surface (Heatmap)
Visualizing the RMSE valley to find the optimal $Q$ matrix parameters. The deep blue region represents the minimum error configuration.
![Optimization Surface](img/heatmap.png)
*(Note: Run `plots/heatmap.gp` to generate this interactive 3D view)*

### 2. Real-Time Tracking & Error Distribution
The filter achieves a perfect Gaussian error distribution centered at zero, demonstrating an unbiased estimator.
![Tracking Dashboard](img/dashboard.png)
*(Note: Run `plots/dashboard.gp` to generate this view)*

## 🛠️ Mathematical Model

The system models a Continuous White Noise Acceleration (CWNA) / Constant Velocity model:

$$x_{k+1} = \begin{bmatrix} 1 & \Delta t \\ 0 & 1 \end{bmatrix} x_k + w_k$$

Where the state vector $x$ contains position and velocity. The core challenge addressed in this repo is tuning the **Process Noise Covariance Matrix ($Q$)**:

$$Q = \begin{bmatrix} \sigma_{pos}^2 & 0 \\ 0 & \sigma_{vel}^2 \end{bmatrix}$$

We iterate over orders of magnitude for $\sigma^2$ to find the global minimum on the error surface.

## 🚀 Getting Started

### Prerequisites
* **GFortran** (GNU Fortran Compiler)
* **Gnuplot** (For rendering the graphs)

### Installation & Execution

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/TU_USUARIO/filtro-kalman-fortran.git](https://github.com/TU_USUARIO/filtro-kalman-fortran.git)
    cd filtro-kalman-fortran
    ```

2.  **Compile the optimizer:**
    ```bash
    gfortran src/kalman_main.f90 -o kalman_solver -O3
    ```

3.  **Run the simulation:**
    ```bash
    ./kalman_solver
    ```
    *This will generate `kalman_continuous.txt` and `error_surface.txt`.*

4.  **Visualize results:**
    ```bash
    gnuplot -p plots/dashboard.gp
    gnuplot -p plots/heatmap.gp
    ```

## 📂 Project Structure

```text
.
├── src/
│   ├── kalman_main.f90      # Main filter logic & optimization loop
│   └── (modules...)
├── plots/
│   ├── dashboard.gp         # Signal tracking & Histogram plot
│   └── heatmap.gp           # 3D Surface optimization plot
├── img/                     # Screenshots for this README
├── README.md
└── .gitignore
