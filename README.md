# rstarGlobal

Replication files for *Global Trends in Interest Rates* by Marco Del Negro, Domenico Giannone, Marc Giannoni, and Andrea Tambalotti.

## Updated r* estimates
The [Excel file](update/Rstar_Vintages.xlsx) contains original and updated estimates (forthcoming) of the U.S. and world real interest rates. These are the estimates displayed in figure 1 of the paper. The below figure shows the latest values.


![Figure 1](update/fig1-Model1_Rshortbar-us.png )
*Fig.1. Trends in Global and U.S. Real Rates: 1870-2016, Baseline Model. Note: The figure plots the posterior median of the trend in the world real interest rate (dashed line) together with its 68 and 95 percent posterior coverage intervals, as well as the posterior median of the trend in the U.S. real rate (dotted line).*

## Required software

These scripts were produced using MATLAB R2017b.

## Installing this repository

Git users are welcome to fork this repository or clone it for local use. Non-Git users will probably find it easiest to download the zip file by clicking on the green `Clone or download` button on the right hand side of this screen, and then clicking "Download ZIP."

## Directory structure
- `figures/`: Figures from the paper's main body.
- `indata/`: Input data.
	- `Data_MY.xlsx`: Middle-aged and young individuals dataset.
	- `DataInflShortLongUpdated.xlsx`: Dataset for baseline model (inflation, short-term rate, long-term rate).
- `results/`: Model results in '.mat' format.
	- `OutputModelX.mat`: See "scripts/" for the naming convention.
- `scripts/`: Scripts for estimation and producing figures
	- `MainModelX_MakeFigures.m`: The "MakeFigures" suffix makes figures for the respective model.
	- `MainModelvar01.mat`: Estimates the baseline model with prior to variance of innovation to the trend component equal to 1.
	- `MainModel1.mat`: Estimates the baseline model.
	- `MainModel1_ReR`: Estimates baseline model under the real exchange rate specification.
	- `MainModel1_unrestr`: Estimates baseline model with unrestricted loading on real rate.
	- `MainModel1_1950.m`: Estimates baseline model starting in 1950.
	- `MainModel1_df50.m`: Estimates baseline model with 50 degrees of freedom.
	- `MainModel1_varX.m`: Estimates baseline model with variance of innovation to trend equal to 1/X.
	- `MainModel2.m`: Estimates convenience yield model.
	- `MainModel3.m`: Estimates consumption model.
- `update/`: Updated versions of figure 1 along with its underlying data.

## How to run the code

Each model is generated using a script titled `MainModelX.m`. Scripts with an underscored suffix (e.g. `MainModel1_1950.m`) correspond to minor modifications to the baseline model; see the appendix for details. 

`estimateAll.m` is the main script for running all the estimation routines. Set `estimateAppendices = 1` to run all specifications. MCMC results are stored as `.mat` files in the `results` folder. Depending on the machine, running the main body specifications for 100,000 draws takes around 20 hours. Reduce the number of draws for faster results.

Paper figures are created using the files with the `_MakeFigures.m` suffix. 

## Data sources

Data on short-term rates, long-term rates, and consumer prices come from the Jordà-Schularick-Taylor Macrohistory Database. Data for Moody's Baa corporate bond yield are available from FRED. Demographic data come from the United States Census Bureau and from the UN World Population Statistics database. Data on corporate spreads come from Gilchrist and Mojon (2018). Please see the paper for additional details on the data sources.
