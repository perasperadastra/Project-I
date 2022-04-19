## Input parameters <a name = "parameters"></a>

The dimension parameters can't be changed in this version.

#### Initial structure

Can be generated three type of lattices: 
  * Simple cubic (structure=1) 
  * Face cubic centered (structure=2)
  * Diamond (structure=3) 
  * Also it will be implemented a read from file subroutine in the next version.

#### Density
The density have units of kg/m^3 . Must specify the number of atoms per side. 


#### Temperature
Temperature is in kelvin units.

#### Initial Velocities

The initial velocities of the particles can be setted up to 0 or to a bimodal distribution (using the temperature of the system). To choose use the **vel_opt** parameter (= 1, bimodal) (= 0, zero initial velocities).


#### Integrator
The integrated algorithms can be the Velocity Verlet with thermostat (= 1) or without thermostat (= 0). 

#### Forces and potential

At last, you must choose the parameters for the force-field ([Lennard-Jones](https://es.wikipedia.org/wiki/Potencial_de_Lennard-Jones) type). 

  * Distance to the zero potential point in hte potential (**σ**): Units in angstroms.
  * Depth of the potential well (**ε**): Units in Kelvins.
  * 

#### Time

To control de time step change the **dt** parameter (in ps units). 
The number of steps that will be computed is selected in the **ntimes** parameter. 
At **tmelt** parameter you can decide some initial steps that initialize the system, represents the time steps of melting the initial structure.

The output information frequency is controlled by the **everyt** parameter. AVOID HIGH RATES OF PRINTIN!!! Printing is a limiting stage in this software.

#### Seed

For reproducibility you can decide a seed for the random number generator that are inside of the software. Remember to be consistent with the seed you use.
