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
//   dt is clamped internally to [1e-6, 0.1] s; non-positive or NaN dt is a no-op.
void bwr_step(BwrHandle handle, double dt);

// ── Operator controls ──────────────────────────────────────────────────────────

// Control rod bank position (BWR bottom-entry convention).
//   0.0 = fully withdrawn (maximum reactivity)
//   1.0 = fully inserted  (minimum reactivity / shutdown)
// Broadcasts the requested fraction to all 137 control blades. For an
// instrumented scram with the spec 2.3 Figure 2.3-13 time-vs-pressure
// curve, prefer bwr_scram() — it latches the scram condition and the
// blades drive to fully inserted following the typical-drive envelope
// (~2.0–2.3 s at 1000 psig, ~3.5 s at 0 psig).
void bwr_set_control_rod_position(BwrHandle handle, double insertion_fraction);

// Latch a reactor protection system scram. All 137 control blades begin
// insertion following the "typical drive" curve from spec 2.3 Figure
// 2.3-13. Idempotent — calling while already scrammed is a no-op.
void bwr_scram(BwrHandle handle);

// Clear the latched scram. Does not retract blades; the operator must
// withdraw them via bwr_set_control_rod_position.
void bwr_scram_reset(BwrHandle handle);

// Core recirculation flow [kg/m²·s]. Legacy interface — routes to the
// recirc pump-speed setpoint: kg_m2_s / rated_mass_flux × 100 % is
// applied to both pumps. The actual mass flux the kernel sees may
// differ during transients due to NPSH protection and pump dynamics;
// for direct pump control use bwr_set_recirc_pump_speed_pct.
void bwr_set_mass_flux(BwrHandle handle, double kg_m2_s);

// Set both recirc pump speed demands [% of rated, 30–102 %]. The pumps
// slew toward this setpoint at the configured speed-slew rate. Below
// 30 % is the turndown floor (spec 2.4 p.3) and collapses speed to zero.
void bwr_set_recirc_pump_speed_pct(BwrHandle handle, double pct);

// Trip a specific recirc pump (idx ∈ {1, 2}). Pump coasts to zero;
// the unaffected loop keeps running. Both pumps off drops core flow to
// the natural-circulation floor (~25 % rated).
void bwr_trip_recirc_pump(BwrHandle handle, int idx);

// Reset all tripped recirc pumps and clear the EOC-RPT latch. Operator
// must still re-establish the pump-speed demand to spin them back up.
void bwr_reset_recirc_pumps(BwrHandle handle);

// Feedwater flow [% of rated steam flow]. Clamped to [0, 200].
// Reactor water level holds at rated power when this is 100 %; below drains,
// above floods. Bypassing or tripping feedwater (set to 0) will drain the
// level steadily and eventually trip the low-level SCRAM.
void bwr_set_feedwater_flow_pct(BwrHandle handle, double pct);

// Current feedwater flow [% of rated steam flow] — the actual delivered
// flow, not the operator demand. Limited by pump-chain capacity:
// one RFP available caps it at 67 %; both RFPs tripped → 0 %.
double bwr_get_feedwater_flow_pct(BwrHandle handle);

// Current heater-train outlet (vessel sparger inlet) temperature [K].
// At rated flow this approaches `rated_inlet_temperature_K` (488.7 K /
// 420 °F by default — Browns Ferry plant rating, spec 2.1); at zero flow
// it drifts to the hotwell temperature (~319 K / 115 °F).
double bwr_get_feedwater_temp_K(BwrHandle handle);

// Count of RFPs (Reactor Feed Pumps) currently running (0..2). Each
// pump is rated 67 % of system capacity (spec 2.6 p.11).
int bwr_get_rfp_running_count(BwrHandle handle);

// Trip a specific RFP. idx is 1-based (1 or 2). Idempotent — re-tripping
// is a no-op; out-of-range indices are ignored. The booster + condensate
// chains continue running; only the specified RFP latches tripped.
void bwr_trip_feedwater_pump(BwrHandle handle, int idx);

// Reset all tripped feedwater pumps (condensate, booster, RFP) and clear
// the staggered suction-trip timers. Used to recover from a transient
// loss-of-FW; operator must still re-establish FW demand via
// bwr_set_feedwater_flow_pct.
void bwr_reset_feedwater_pumps(BwrHandle handle);

// Turbine steam admission valve [0 – 100 %]. Commands the TCV / TSV
// demand; the actual valve position slews at ~100 %/s toward this.
// A turbine trip (bwr_trip_turbine) overrides with a 0.1 s fast-close
// ramp regardless of demand.
void bwr_set_turbine_valve(BwrHandle handle, double pct);

// Steam bypass valve demand [0 – 100 %]. The valve stroke is the value
// reported by bwr_get_bypass_valve_pct; the effective steam relief is
// capped at 25 % of rated steam flow even at 100 % stroke
// (spec 2.5 p.11).
void bwr_set_bypass_valve(BwrHandle handle, double pct);

// Latch a TSV "fast closure" turbine trip. The four turbine stop valves
// ramp to fully closed over 100 ms (spec 2.5 p.12); the turbine_tripped
// flag latches once they reach zero. Idempotent.
void bwr_trip_turbine(BwrHandle handle);

// Clear the latched turbine trip. Operator must re-open the TCVs via
// bwr_set_turbine_valve afterwards.
void bwr_reset_turbine_trip(BwrHandle handle);

// Latch the NSSSS MSIV isolation signal. All four MSIVs stroke shut
// over ~4 s ([BWR/4 typical 3-5 s]; spec 2.5 p.20). Idempotent.
void bwr_close_msivs(BwrHandle handle);

// Clear the NSSSS MSIV isolation signal. MSIVs stroke back open over
// the same stroke time. Should only be called once the underlying
// isolation condition has cleared.
void bwr_open_msivs(BwrHandle handle);

// Adjust the sustained operator reactivity perturbation [pcm].
//   Positive = reactivity addition (power up).
//   Negative = reactivity removal  (power down).
//   Each call ADDS to the running total. The perturbation is reapplied each
//   step by scaling nu*sigma_f, so the effect persists until cleared with
//   bwr_reset_reactivity().
void bwr_apply_reactivity(BwrHandle handle, double rho_pcm);

// Zero the sustained operator reactivity perturbation.
void bwr_reset_reactivity(BwrHandle handle);

// Current accumulated operator reactivity perturbation [pcm].
double bwr_get_applied_reactivity_pcm(BwrHandle handle);

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

// Water level above the top of active fuel [m]. Negative values indicate
// fuel uncovery (water level has dropped below the top of the active core).
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

// Turbine rotor speed [rpm]. BWR/4 tandem-compound 4-flow turbine drives
// a 4-pole synchronous generator: rated speed = 1800 rpm @ 60 Hz.
double bwr_get_turbine_speed_rpm(BwrHandle handle);

// Mechanical shaft power [W]. Equals (eta_turbine × pressure_factor ×
// tcv_flow_norm × power_current_W). Goes to zero when the turbine is
// tripped or no steam is flowing through the TCVs.
double bwr_get_turbine_mech_W(BwrHandle handle);

// Steam fraction through the turbine control valves (drives the rotor).
//   1.0 = rated steam flow at full TCV / rated pressure.
double bwr_get_tcv_flow_norm(BwrHandle handle);

// Steam fraction through the bypass valves (dumps to condenser, no shaft
// work). Hardware-capped at 0.25 of rated even at 100 % stroke.
double bwr_get_bpv_flow_norm(BwrHandle handle);

// Indicated grid frequency [Hz], derived from turbine_speed_rpm / rated × 60.
// At rated 1800 rpm this reads 60.0 Hz. Useful for the sync window UI.
double bwr_get_grid_frequency_Hz(BwrHandle handle);

// ── Manual generator-breaker controls ─────────────────────────────────────────
// The operator manually closes the generator breaker when the turbine is
// near rated speed (±sync_window_rpm = 15 RPM around 1800). Outside that
// window the call silently no-ops; poll bwr_get_sync_ready() to confirm
// the breaker may close.
void bwr_sync_generator(BwrHandle handle);

// Open the generator breaker (load reject / desync). Drops the turbine
// back to free-spinning; without operator action on the TCV the rotor
// will overspeed up to free_speed_overshoot × rated_rpm.
void bwr_unsync_generator(BwrHandle handle);

// 1 if the generator breaker is closed and rotor is locked to grid
// frequency. 0 otherwise. Goes to 0 on turbine trip and on reverse-
// power protection (mechanical power < 50 % of bearing loss).
int bwr_get_generator_synced(BwrHandle handle);

// 1 when the breaker may be safely closed: not tripped, not already
// synced, and rotor speed within ±sync_window_rpm of rated.
int bwr_get_sync_ready(BwrHandle handle);

// Current TCV / TSV stroke position [0 – 100 %]. Differs from the
// commanded demand during the 100 ms TSV fast-close ramp on trip.
double bwr_get_turbine_valve_pct(BwrHandle handle);

// Current bypass valve stroke position [0 – 100 %].
double bwr_get_bypass_valve_pct(BwrHandle handle);

// Average MSIV open fraction across the four steam lines [0 – 1].
// 1.0 = all MSIVs fully open; 0.0 = all four MSLs fully isolated.
double bwr_get_msiv_open_frac(BwrHandle handle);

// Number of SRVs currently lifted (0 – 11).
int bwr_get_srv_count_open(BwrHandle handle);

// Total SRV blowdown flow [kg/s]. Each lifted SRV contributes its
// per-valve capacity (~106 kg/s, typical BWR/4).
double bwr_get_srv_flow_kg_s(BwrHandle handle);

// Current core mass flux fed to the two-phase model [kg/m²·s]. At
// rated conditions this is the value passed to bwr_reactor_create.
double bwr_get_mass_flux(BwrHandle handle);

// Total core mass flow [kg/s], from the jet-pump M-ratio model:
// core_flow = drive_flow × (1 + M). Floors at the natural-circulation
// fraction (~25 % rated) when both pumps are tripped.
double bwr_get_core_flow_kg_s(BwrHandle handle);

// Core flow as a fraction of rated [%]. 100 % at rated, ~25 % on
// natural circulation, intermediate values during pump runback.
double bwr_get_core_flow_pct(BwrHandle handle);

// Current speed of recirc pump idx [% of rated, 1-based]. Returns −1 if
// idx is out of range.
double bwr_get_recirc_pump_speed_pct(BwrHandle handle, int idx);

// Current control rod bank insertion fraction [0.0 – 1.0].
// Returns the average insertion fraction across all 137 blades. For
// uniform-bank operation this is the same as any individual blade's
// position.
double bwr_get_control_rod_position(BwrHandle handle);

// Volume-averaged fuel burnup [MWd/kgU].
double bwr_get_avg_burnup(BwrHandle handle);

// Volume-averaged coolant temperature [K]. The bulk subcooled temperature
// the convection helper uses to compute fuel-to-coolant heat transfer.
double bwr_get_avg_coolant_temp_K(BwrHandle handle);

// Total steam flow leaving the dome, normalised to rated (0..1+).
//   = bwr_get_tcv_flow_norm + bwr_get_bpv_flow_norm
double bwr_get_steam_flow_norm(BwrHandle handle);

// Total steam flow leaving the dome [kg/s]. At rated this is ~1638 kg/s.
double bwr_get_steam_flow_kg_s(BwrHandle handle);

// Feedwater mass flow into the vessel [kg/s]. Pair with the steam flow
// to monitor the inventory balance that drives reactor level.
double bwr_get_feedwater_flow_kg_s(BwrHandle handle);

// Simulated time elapsed since reactor_create [s].
double bwr_get_time_s(BwrHandle handle);

// Number of bwr_step calls since reactor_create.
int bwr_get_step_count(BwrHandle handle);

// Rated thermal power [W] (the value passed to bwr_reactor_create). Use
// this as the denominator for relative-power displays so the scale stays
// independent of the operating point.
double bwr_get_power_rated_W(BwrHandle handle);

// Turbine trip latch (1 = TSV fast-close latched, 0 = normal).
int bwr_get_turbine_tripped(BwrHandle handle);

// SCRAM latch (1 = rod drives commanded to full insertion, 0 = normal).
int bwr_get_scram_latched(BwrHandle handle);

// MSIV close latch (1 = NSSSS isolation commanded / in progress, 0 = open).
int bwr_get_msiv_close_latched(BwrHandle handle);

// Rod Block Monitor latch (1 = withdrawal frozen due to short reactor
// period, 0 = motion permitted). The block clears automatically once
// period recovers above the configured setpoint (default 30 s). Insertion
// is never blocked.
int bwr_get_rod_block_active(BwrHandle handle);

// Operator rod-demand fraction [0..1]. Updates instantly when the operator
// commands a position; the actual blade insertion (bwr_get_control_rod_position)
// slews toward this at the hardware speed. When the two differ, rods are
// either still moving or held by the RBM.
double bwr_get_rod_demand(BwrHandle handle);

// ── Feedwater control mode ────────────────────────────────────────────────────
// 3-element controller toggle. When ON, demand_flow_pct is overwritten
// each tick to match ṁ_steam + a level-error trim, so feedwater tracks
// steam outflow automatically. When OFF, the operator demand from
// bwr_set_feedwater_flow_pct is used verbatim. The controller also
// auto-engages once dome pressure crosses the boiling-active threshold
// (~2 MPa) unless explicitly disabled here.
void bwr_set_fw_controller_enabled(BwrHandle handle, int on);
int  bwr_get_fw_controller_enabled(BwrHandle handle);

// ── Residual Heat Removal (RHR) ───────────────────────────────────────────────
// 4-loop RHR system (spec 10.4). Modes:
//   0 STANDBY           — pumps idle, valves lined up for LPCI auto-init
//   1 SHUTDOWN_COOLING  — coolant direct-cooling via RHR HX (permissive
//                         dome P < 135 psig)
//   2 SUPP_POOL_COOLING — suppression pool → HX → pool
//   3 CONTAINMENT_SPRAY — pool → HX → drywell/wetwell spray (permissive
//                         drywell P > 1.7 psig)
//   4 LPCI              — pool → recirc loop discharge (ECCS injection)
// Permissive checks reject mode changes silently; poll the getter to
// confirm the mode actually changed.
#define BWR_RHR_MODE_STANDBY           0
#define BWR_RHR_MODE_SHUTDOWN_COOLING  1
#define BWR_RHR_MODE_SUPP_POOL_COOLING 2
#define BWR_RHR_MODE_CONTAINMENT_SPRAY 3
#define BWR_RHR_MODE_LPCI              4

void bwr_set_rhr_mode(BwrHandle handle, int mode);
int  bwr_get_rhr_mode(BwrHandle handle);

// Start (on != 0) or stop (on == 0) RHR loop pump idx (1..4).
void bwr_set_rhr_pump(BwrHandle handle, int idx, int on);

// Count of RHR pumps currently running (not tripped, command = on).
int bwr_get_rhr_pumps_running(BwrHandle handle);

// Total heat removed by RHR this tick [W]. Zero in STANDBY. Goes into
// the reactor coolant in SDC, into the suppression pool in the other
// running modes.
double bwr_get_rhr_total_heat_W(BwrHandle handle);

// RHR heat-exchanger outlet temperature [K]. Useful diagnostic when
// pool-cooling: pool T should converge to (T_in − heat / mc_dot).
double bwr_get_rhr_hx_outlet_T_K(BwrHandle handle);

// Suppression pool bulk temperature [K]. Heats up from SRV blowdown,
// cools via RHR pool-cooling mode.
double bwr_get_supp_pool_T_K(BwrHandle handle);

// ── Condensate loop (water-side, post-hotwell) ───────────────────────────────
// Hotwell water level [m]. Mid-band normal ≈ 2.0 m; low-level alarm
// triggers makeup demand from CST. Tracks the integrated steam vs FW
// mass balance through the condenser.
double bwr_get_hotwell_level_m(BwrHandle handle);

// Condensate Storage Tank water level [m]. ~9 m at startup full; 2.45 m
// floor protects the 150 kgal HPCI/RCIC reserve.
double bwr_get_cst_level_m(BwrHandle handle);

// Hotwell makeup flow demand [kg/s]. Non-zero when hotwell level is
// below the low-level threshold; supplied from the CST via the
// transfer pumps.
double bwr_get_cond_makeup_kg_s(BwrHandle handle);

// Dissolved-O₂ in the condensate [ppb]. Normal vacuum: ~8 ppb (decays
// toward normal under SJAE removal). Lost vacuum: rises rapidly toward
// the high-alarm threshold (50 ppb).
double bwr_get_cond_o2_ppb(BwrHandle handle);

// Print a multi-line state summary to stdout. Useful when the frontend
// can pipe the sim's stdout into its console.
void bwr_print_state_summary(BwrHandle handle);

// ── Safety trips ───────────────────────────────────────────────────────────────
// bwr_get_trip_flags returns a bitmask of currently-tripped safety thresholds.
// A non-zero return means the reactor satisfies at least one SCRAM condition;
// the GDExtension wrapper latches the SCRAM and emits the signal. Callers may
// also poll between steps to drive pre-trip alarm indicators.
#define BWR_TRIP_FUEL_TEMP_HIGH    (1 << 0)  // peak fuel T > 1400 °C
#define BWR_TRIP_POWER_HIGH        (1 << 1)  // total power > 130 % rated
#define BWR_TRIP_SHORT_PERIOD      (1 << 2)  // 0 < period < 5 s
#define BWR_TRIP_LEVEL_LOW         (1 << 3)  // level < +0.3 m above top of fuel (L3 scram)
#define BWR_TRIP_PRESSURE_HIGH     (1 << 4)  // pressure > 7.45 MPa (~1080 psig high-P scram)
#define BWR_TRIP_PRESSURE_LOW      (1 << 5)  // pressure < 5.50 MPa (~798 psig low-MSL iso)
#define BWR_TRIP_CHFR_LOW          (1 << 6)  // critical heat flux ratio < 1

int bwr_get_trip_flags(BwrHandle handle);

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

// Debug: returns center-cell thermal (group-2) absorption XS and the
// internal mg_state k_eff. Used to verify XS priming and Chebyshev fix.
void bwr_debug_xs(BwrHandle handle, double *sigma_a2, double *mg_keff);

#ifdef __cplusplus
}
#endif

#endif // NUCLEAR_PHYSICS_BWR_H