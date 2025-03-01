# pyLBL

## Installation
Make sure that you have the most recent version of `pip`, then run
the following command in the base directory of the repository:

```python
pip install .
```

## High-level API
The calculation of absorption coefficients requires a `Spectroscopy` object.  If a
`database` parameter is provided upon initialization, the object will create a local
SQLite3 database containing the data necessary to calculate the absorption for the
list of input molecules if a file by that name does not exist.  If the file does exist,
the object will assume that file is a SQLite3 database that conforms to the schema detailed
below.  After intialization, the specral parameters can be loaded into memory by calling
the `load_spectral_inputs` method.

```python
from pyLBL import Spectroscopy

spectroscopy = Spectroscopy(["H2O", "CO2"], database="local_database.sqlite")
spectroscopy.load_spectral_inputs()
spectroscopy.list_molecules()
```


### Spectral database schema
Spectral data is downloaded from the web and stored as tables in a local SQLite3 database.
The data for each component of the absorption is stored as follows:


#### HITRAN line parameter tables
Line parameters for each molecule are taken from [HITRAN](https://hitran.org) and stored in
tables with columns:
| CH4_lines       |
|-----------------|
| id INTEGER      |
| iso INTEGER     |
| v REAL          |
| s REAL          |
| en REAL         |
| d_air REAL      |
| gamma_air REAL  |
| gamma_self REAL |
| n_air REAL      |
| n_self REAL     |


#### Total Internal Partition Sums
Total internal parition sum parameters for each molecule are taken from
[https://doi.org/10.1016/j.jqsrt.2017.03.045](https://doi.org/10.1016/j.jqsrt.2017.03.045)
and stored in tables with columns:

| CH4_tips         |
|------------------|
| temperature REAL |
| Q_1 REAL         |
| Q_2 REAL         |
| ...              |


#### Absorption cross sections
Absorption cross-sections for each molecule are taken from [HITRAN](https://hitran.org)
and stored in a set of related tables with columns:

| N2O_cross_section_band_parameters |
|-----------------------------------|
| band_id INTEGER PRIMARY KEY       |
| lower_bound REAL                  |
| upper_bound REAL                  |
| size INTEGER                      |

| N2O_cross_section_temperatures     |
|------------------------------------|
| temperature_id INTEGER PRIMARY KEY |
| temperature REAL UNIQUE            |

| N2O_cross_section_pressures     |
|---------------------------------|
| pressure_id INTEGER PRIMARY KEY |
| pressure REAL UNIQUE            |

| N2O_cross_sections                                                            |
|-------------------------------------------------------------------------------|
| wavenumber REAL                                                               |
| cross_section REAL                                                            |
| band INTEGER REFERENCES N2O_cross_section_band_parameters(band_id)            |
| temperature INTEGER REFERENCES N2O_cross_section_temperatures(temperature_id) |
| pressure INTEGER REFERENCES N2O_cross_section_pressures(pressure_id)          |


#### Collision-induced Absorption.
Collision-induced absorption parameters are also taken from [HITRAN](https://hitran.org)
and stored a in a set of related tables with columns:

| CO2_CH4_cia_band_parameters |
|-----------------------------|
| band_id INTEGER PRIMARY KEY |
| lower_bound REAL            |
| upper_bound REAL            |
| size INTEGER                |

| CO2_CH4_cia_temperatures           |
|------------------------------------|
| temperature_id INTEGER PRIMARY KEY |
| temperature REAL UNIQUE            |

| CO2_CH4_cia_cross_sections                                              |
|-------------------------------------------------------------------------|
| wavenumber REAL                                                         |
| cross_section REAL                                                      |
| band INTEGER REFERENCES CO2_CH4_cia_band_parameters(band_id)            |
| temperature INTEGER REFERENCES CO2_CH4_cia_temperatures(temperature_id) |


#### Water vapor continuum absorption
Water-vapor continuum absorption coefficient parameters are stored in tables with columns:

| H2O_foreign_continuum |
|-----------------------|
| wavenumber REAL       |
| c REAL                |
| t REAL                |

| H2O_self_continuum |
|--------------------|
| wavenumber REAL    |
| c REAL             |
| t REAL             |


#### Ozone continuum.
Ozone continuum absorption coefficient parameters are stored in a table with columns:

| O3_continuum       |
|--------------------|
| wavenumber REAL    |
| cross_section REAL |
