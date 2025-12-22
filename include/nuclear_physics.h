// include/nuclear_physics.h
//
// C interface for the nuclear physics library
//
#ifndef NUCLEAR_PHYSICS_H
#define NUCLEAR_PHYSICS_H

#ifdef __cplusplus
extern "C" {
#endif

// Opaque reactor handle
typedef void* ReactorHandle;

// Create reactor simulation
// Parameters:
//   nx, ny, nz: Grid dimensions
//   dx, dy, dz: Grid spacing [m]
// Returns: Handle to reactor instance
ReactorHandle reactor_create(int nx, int ny, int nz, 
                             double dx, double dy, double dz);

// Destroy reactor simulation
void reactor_destroy(ReactorHandle reactor);

// Set fuel region with nuclear properties
// Parameters:
//   i1, i2, j1, j2, k1, k2: Region bounds (inclusive)
//   power_density: Nominal power density [W/m³]
//   enrichment: U-235 enrichment fraction [0-1]
void reactor_set_fuel_region(ReactorHandle reactor,
                             int i1, int i2, int j1, int j2, int k1, int k2,
                             double power_density, double enrichment);

// Set coolant (water) region
// Parameters:
//   inlet_temp: Inlet temperature [K]
//   mass_flow: Mass flow rate [kg/s]
void reactor_set_coolant_region(ReactorHandle reactor,
                                int i1, int i2, int j1, int j2, int k1, int k2,
                                double inlet_temp, double mass_flow);

// Set control rod position
// Parameters:
//   insertion_fraction: 0 = fully withdrawn, 1 = fully inserted
void reactor_set_control_rods(ReactorHandle reactor, double insertion_fraction);

// Perform one time step
// Parameters:
//   dt: Time step size [s]
void reactor_step(ReactorHandle reactor, double dt);

// Get temperature field
// Parameters:
//   T_out: Output array (must be allocated with size nx*ny*nz)
void reactor_get_temperature(ReactorHandle reactor, double* T_out);

// Get pressure field
// Parameters:
//   p_out: Output array (must be allocated with size nx*ny*nz)
void reactor_get_pressure(ReactorHandle reactor, double* p_out);

// Get power density field
// Parameters:
//   power_out: Output array (must be allocated with size nx*ny*nz)
void reactor_get_power(ReactorHandle reactor, double* power_out);

// Get total reactor thermal power [W]
double reactor_get_total_power(ReactorHandle reactor);

// Get velocity field components
// Parameters:
//   vx_out, vy_out, vz_out: Output arrays (each size nx*ny*nz)
void reactor_get_velocity(ReactorHandle reactor, 
                         double* vx_out, double* vy_out, double* vz_out);

// Get maximum stable time step [s]
double reactor_get_max_dt(ReactorHandle reactor);

#ifdef __cplusplus
}
#endif

#endif // NUCLEAR_PHYSICS_H