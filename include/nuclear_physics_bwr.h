// include/nuclear_physics_bwr.h
//
// C interface for the full BWR simulation.
// Wraps coupled multigroup-diffusion / heat-transfer / two-phase / burnup physics.
//
// Quick start:
//   BwrHandle bwr = bwr_reactor_create(
//       20, 20, 20,                  // nx, ny, nz
//       0.2375, 0.2375, 0.1855,      // dx, dy, dz [m]  (4.75m / 20, 3.71m / 20)
//       2381e6,                      // rated thermal power [W]
//       551.0,                       // inlet temperature [K]
//       7.14e6,                      // operating pressure [Pa]
//       1500.0,                      // core mass flux [kg/m²·s]
//       0.035,                       // U-235 enrichment
//       3.71,                        // active core height [m]
//       4.75                         // core diameter [m]
//   );
//   bwr_solve_steady_state(bwr);     // converge initial conditions (~50 iters)
//   bwr_step(bwr, 0.02);             // advance 20 ms (50 Hz sim tick)
//   double keff = bwr_get_keff(bwr);
//   bwr_reactor_destroy(bwr);
//
// Field memory layout (bwr_get_*_field functions):
//   Arrays are Fortran column-major: i varies fastest.
//   Element at 1-based cell (i, j, k) is out[i-1 + nx*(j-1 + ny*(k-1))].
//   Allocate nx*ny*nz doubles before calling any field getter.
//
#ifndef NUCLEAR_PHYSICS_BWR_H
#define NUCLEAR_PHYSICS_BWR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void* BwrHandle;

// ── Lifecycle ──────────────────────────────────────────────────────────────────
BwrHandle bwr_reactor_create(int nx, int ny, int nz,
                              double dx, double dy, double dz,
                              double power_rated_W,
                              double inlet_temp_K,
                              double pressure_Pa,
                              double mass_flux,
                              double enrichment,
                              double core_height_m,
                              double core_diameter_m);

void bwr_reactor_destroy(BwrHandle handle);

// ── Simulation control ─────────────────────────────────────────────────────────
// bwr_solve_steady_state: run up to 50 coupled outer iterations to convergence.
//   Call once after create, before the transient loop.
void bwr_solve_steady_state(BwrHandle handle);

// bwr_step: advance the coupled simulation by dt seconds.
//   Recommended dt <= 0.02 s (50 Hz). Each call runs:
//     neutronics → burnup → heat transfer → two-phase → XS feedback →
//     instrumentation → pressure dynamics → turbine
void bwr_step(BwrHandle handle, double dt);

// ── Operator controls ──────────────────────────────────────────────────────────

// Control rod bank position (BWR bottom-entry convention).
//   0.0 = fully withdrawn (maximum reactivity)
//   1.0 = fully inserted  (minimum reactivity / shutdown)
void bwr_set_control_rod_position(BwrHandle handle, double insertion_fraction);

// Core coolant mass flux [kg/m²·s]. Clamped >= 100.
void bwr_set_mass_flux(BwrHandle handle, double kg_m2_s);

// Turbine steam admission valve [0 – 100 %].
void bwr_set_turbine_valve(BwrHandle handle, double pct);

// Steam bypass valve [0 – 100 %].
void bwr_set_bypass_valve(BwrHandle handle, double pct);

// Instantaneous reactivity insertion [pcm].
//   Positive = reactivity addition (power up).
//   Negative = reactivity removal  (power down).
//   Implemented by scaling nu*sigma_f; accumulates on repeated calls.
void bwr_apply_reactivity(BwrHandle handle, double rho_pcm);

// ── Scalar instrument readings ─────────────────────────────────────────────────

// Effective neutron multiplication factor.
double bwr_get_keff(BwrHandle handle);

// Core reactivity [(k_eff - 1) / k_eff * 1e5] in pcm.
double bwr_get_reactivity_pcm(BwrHandle handle);

// Total thermal power [W].
double bwr_get_total_power_W(BwrHandle handle);

// Peak fuel temperature [K].
double bwr_get_max_fuel_temp_K(BwrHandle handle);

// Volume-averaged void fraction [%] (0 – 100).
double bwr_get_avg_void_fraction(BwrHandle handle);

// Minimum critical heat flux ratio in the core. Values < 1.0 indicate CHF.
double bwr_get_min_chfr(BwrHandle handle);

// Current reactor operating pressure [Pa].
double bwr_get_pressure_Pa(BwrHandle handle);

// Water level above core floor [m]. Negative values indicate uncovery.
double bwr_get_reactor_level_m(BwrHandle handle);

// Saturation temperature at current pressure [K].
double bwr_get_sat_temperature_K(BwrHandle handle);

// Average Power Range Monitor [% of rated power]. Clamped to [0, 200].
double bwr_get_aprm_pct(BwrHandle handle);

// Reactor period [s]. Large positive values indicate a stable/subcritical reactor.
//   Values in (0, 5) indicate prompt supercriticality risk.
double bwr_get_reactor_period_s(BwrHandle handle);

// Net turbine electrical output [W]. Zero when tripped.
double bwr_get_turbine_power_W(BwrHandle handle);

// Turbine rotor speed [rpm]. Synchronised grid speed is 3600 rpm.
double bwr_get_turbine_speed_rpm(BwrHandle handle);

// Current control rod bank insertion fraction [0.0 – 1.0].
double bwr_get_control_rod_position(BwrHandle handle);

// Volume-averaged fuel burnup [MWd/kgU].
double bwr_get_avg_burnup(BwrHandle handle);

// ── 3-D field outputs ──────────────────────────────────────────────────────────
// Caller must pre-allocate out with nx*ny*nz doubles.
// Memory order: i varies fastest (Fortran column-major).
// Cell (i,j,k) [1-based] -> out[i-1 + nx*(j-1 + ny*(k-1))].

// Power density [W/m³].
void bwr_get_power_density(BwrHandle handle, double* out);

// Fuel/coolant temperature [K].
void bwr_get_temperature(BwrHandle handle, double* out);

// Local void fraction [0 – 1].
void bwr_get_void_fraction(BwrHandle handle, double* out);

// Local fuel burnup [MWd/kgU].
void bwr_get_burnup_field(BwrHandle handle, double* out);

// Xe-135 concentration [atoms/barn-cm].
void bwr_get_xenon(BwrHandle handle, double* out);

// ── Grid dimensions ────────────────────────────────────────────────────────────
int bwr_get_nx(BwrHandle handle);
int bwr_get_ny(BwrHandle handle);
int bwr_get_nz(BwrHandle handle);

#ifdef __cplusplus
}
#endif

#endif // NUCLEAR_PHYSICS_BWR_H
