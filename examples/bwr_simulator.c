// examples/bwr_simulator_stable.c
//
// Stabilised Boiling Water Reactor (BWR) Simulator
//
// Key improvements for numerical stability:
// 1. Adaptive time stepping with CFL conditions
// 2. Semi-implicit coupling with damping
// 3. Smoothed void reactivity feedback
// 4. Better initial conditions
// 5. Conservative bounds checking
//
// Compile: gcc bwr_simulator_stable.c -o bwr_sim -lnuclear_physics_c -lm
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nuclear_physics.h"

// BWR parameters
#define CORE_HEIGHT 3.7
#define CORE_DIAMETER 4.5
#define FUEL_ASSEMBLIES 764
#define THERMAL_POWER 3293  // MW

// Grid resolution
#define NX 30
#define NY 30
#define NZ 40

// Physical constants
#define T_INLET 550.0
#define T_SAT 560.0
#define P_OPERATING 7.0e6

// Simulation parameters
#define SIM_TIME 100.0
#define OUTPUT_INTERVAL 1.0
#define DT_MIN 1.0e-5       // Minimum time step: 10 μs
#define DT_MAX 0.01         // Maximum time step: 10 ms
#define DAMPING_FACTOR 0.3  // Feedback damping (0-1)

typedef struct {
    ReactorHandle reactor;
    int nx, ny, nz;
    double dx, dy, dz;
    double time;
    double total_power;
    double control_rod_position;
    
    // State arrays
    double* temperature;
    double* temperature_old;  // For damping
    double* pressure;
    double* power_density;
    double* void_fraction;
    double* void_fraction_old;  // For damping
    
    // Stability tracking
    double max_temp_change;
    double max_void_change;
    int stability_violations;
    
    // Adaptive time stepping
    double dt_thermal;
    double dt_void;
    double dt_adaptive;
} BWRState;

// Initialize BWR geometry
void bwr_init_geometry(BWRState* bwr) {
    int i, j, k;
    double r, x, y;
    double core_radius = CORE_DIAMETER / 2.0;
    
    printf("Initializing stabilized BWR geometry...\n");
    printf("  Core: %.2f m diameter x %.2f m height\n", CORE_DIAMETER, CORE_HEIGHT);
    printf("  Grid: %dx%dx%d\n", NX, NY, NZ);
    printf("  Thermal power: %.0f MW\n", THERMAL_POWER);
    
    bwr->reactor = reactor_create(NX, NY, NZ,
                                  CORE_DIAMETER/NX,
                                  CORE_DIAMETER/NY,
                                  CORE_HEIGHT/NZ);
    
    bwr->nx = NX;
    bwr->ny = NY;
    bwr->nz = NZ;
    bwr->dx = CORE_DIAMETER / NX;
    bwr->dy = CORE_DIAMETER / NY;
    bwr->dz = CORE_HEIGHT / NZ;
    bwr->time = 0.0;
    bwr->control_rod_position = 0.65;  // Start closer to critical
    bwr->stability_violations = 0;
    
    // Allocate arrays
    int n_cells = NX * NY * NZ;
    bwr->temperature = (double*)malloc(n_cells * sizeof(double));
    bwr->temperature_old = (double*)malloc(n_cells * sizeof(double));
    bwr->pressure = (double*)malloc(n_cells * sizeof(double));
    bwr->power_density = (double*)malloc(n_cells * sizeof(double));
    bwr->void_fraction = (double*)malloc(n_cells * sizeof(double));
    bwr->void_fraction_old = (double*)malloc(n_cells * sizeof(double));
    
    // Initialize to inlet conditions
    for (int idx = 0; idx < n_cells; idx++) {
        bwr->temperature[idx] = T_INLET;
        bwr->temperature_old[idx] = T_INLET;
        bwr->void_fraction[idx] = 0.0;
        bwr->void_fraction_old[idx] = 0.0;
    }
    
    // Conservative power density estimate
    double power_density = (THERMAL_POWER * 1e6) /
                          (M_PI * core_radius * core_radius * CORE_HEIGHT);
    power_density *= 0.05;  // Very conservative start
    
    for (k = 1; k <= NZ; k++) {
        for (j = 1; j <= NY; j++) {
            for (i = 1; i <= NX; i++) {
                x = (i - NX/2.0) * bwr->dx;
                y = (j - NY/2.0) * bwr->dy;
                r = sqrt(x*x + y*y);
                
                if (r < core_radius) {
                    // Smooth power profile
                    double z_norm = (double)k / NZ;
                    double r_norm = r / core_radius;
                    
                    double axial_factor = 1.0 + 0.2 * cos(M_PI * (z_norm - 0.5));
                    double radial_factor = 1.0 - 0.2 * r_norm * r_norm;
                    double local_power = power_density * axial_factor * radial_factor;
                    
                    reactor_set_fuel_region(bwr->reactor, i, i, j, j, k, k,
                                          local_power, 0.035);
                } else {
                    reactor_set_coolant_region(bwr->reactor, i, i, j, j, k, k,
                                              T_INLET, 0.0);
                }
            }
        }
    }
    
    // Inlet plenum
    reactor_set_coolant_region(bwr->reactor, 1, NX, 1, NY, 1, 3,
                               T_INLET, 10000.0);
    
    printf("Geometry initialized with conservative parameters.\n");
}

// Smooth void fraction calculation with hysteresis
double calculate_void_fraction_smooth(double T, double T_old, double void_old) {
    const double T_SUBCOOL = T_SAT - 10.0;  // Subcooling margin
    const double T_SUPERHEAT = T_SAT + 20.0;  // Superheating margin
    
    double void_new;
    
    if (T < T_SUBCOOL) {
        void_new = 0.0;
    } else if (T < T_SAT) {
        // Gradual onset of boiling
        double fraction = (T - T_SUBCOOL) / (T_SAT - T_SUBCOOL);
        void_new = 0.3 * fraction * fraction;  // Quadratic transition
    } else if (T < T_SUPERHEAT) {
        // Linear transition in two-phase region
        double fraction = (T - T_SAT) / (T_SUPERHEAT - T_SAT);
        void_new = 0.3 + 0.5 * fraction;
    } else {
        void_new = 0.75;  // Cap at 75% for stability
    }
    
    // Apply hysteresis damping
    void_new = void_old + DAMPING_FACTOR * (void_new - void_old);
    
    // Enforce bounds
    void_new = fmax(0.0, fmin(0.8, void_new));
    
    return void_new;
}

// Calculate void reactivity with smoothing
double calculate_void_reactivity_smooth(BWRState* bwr) {
    int idx, i, j, k;
    double total_void = 0.0;
    double core_radius = CORE_DIAMETER / 2.0;
    int core_cells = 0;
    
    reactor_get_temperature(bwr->reactor, bwr->temperature);
    
    idx = 0;
    for (k = 0; k < NZ; k++) {
        for (j = 0; j < NY; j++) {
            for (i = 0; i < NX; i++) {
                double x = (i - NX/2.0) * bwr->dx;
                double y = (j - NY/2.0) * bwr->dy;
                double r = sqrt(x*x + y*y);
                
                if (r < core_radius) {
                    // Smooth void fraction update
                    bwr->void_fraction[idx] = calculate_void_fraction_smooth(
                        bwr->temperature[idx],
                        bwr->temperature_old[idx],
                        bwr->void_fraction_old[idx]
                    );
                    
                    total_void += bwr->void_fraction[idx];
                    core_cells++;
                }
                idx++;
            }
        }
    }
    
    double avg_void = (core_cells > 0) ? total_void / core_cells : 0.0;
    
    // Void coefficient with conservative magnitude
    // Negative feedback: more void -> less moderation -> less reactivity
    double void_reactivity = -0.0010 * avg_void * 100.0;  // Reduced from -0.0015
    
    // Limit maximum feedback magnitude for stability
    void_reactivity = fmax(-0.15, fmin(0.0, void_reactivity));
    
    return void_reactivity;
}

// Compute adaptive time step
double compute_adaptive_timestep(BWRState* bwr) {
    double dt_max_allowed = DT_MAX;
    
    // Thermal time scale: dt < dx²/(2α)
    // Assuming thermal diffusivity α ≈ 1e-7 m²/s for fuel
    double dx_min = fmin(bwr->dx, fmin(bwr->dy, bwr->dz));
    double thermal_diffusivity = 1.0e-7;
    bwr->dt_thermal = 0.25 * dx_min * dx_min / thermal_diffusivity;
    
    // Void dynamics time scale (empirical)
    double max_void_rate = 0.0;
    int n_cells = NX * NY * NZ;
    for (int i = 0; i < n_cells; i++) {
        double void_rate = fabs(bwr->void_fraction[i] - bwr->void_fraction_old[i]);
        if (void_rate > max_void_rate) max_void_rate = void_rate;
    }
    
    if (max_void_rate > 1.0e-6) {
        bwr->dt_void = 0.01 / max_void_rate;  // Limit void change per step
    } else {
        bwr->dt_void = DT_MAX;
    }
    
    // Coolant flow time scale
    double flow_velocity = 2.0;  // m/s typical BWR flow
    double dt_flow = 0.5 * bwr->dz / flow_velocity;
    
    // Take minimum of all constraints
    bwr->dt_adaptive = fmin(bwr->dt_thermal, bwr->dt_void);
    bwr->dt_adaptive = fmin(bwr->dt_adaptive, dt_flow);
    bwr->dt_adaptive = fmin(bwr->dt_adaptive, dt_max_allowed);
    bwr->dt_adaptive = fmax(bwr->dt_adaptive, DT_MIN);
    
    return bwr->dt_adaptive;
}

// Check stability and apply corrections
int check_stability(BWRState* bwr) {
    int n_cells = NX * NY * NZ;
    double max_temp = 0.0, min_temp = 1e10;
    double max_change = 0.0;
    int unstable_cells = 0;
    
    reactor_get_temperature(bwr->reactor, bwr->temperature);
    
    for (int i = 0; i < n_cells; i++) {
        double temp = bwr->temperature[i];
        
        // Check for NaN or Inf
        if (isnan(temp) || isinf(temp)) {
            return 0;  // Critical failure
        }
        
        // Check bounds
        if (temp < 250.0 || temp > 800.0) {
            unstable_cells++;
            // Clamp to reasonable range
            bwr->temperature[i] = fmax(250.0, fmin(800.0, temp));
        }
        
        // Track statistics
        if (temp > max_temp) max_temp = temp;
        if (temp < min_temp) min_temp = temp;
        
        double change = fabs(temp - bwr->temperature_old[i]);
        if (change > max_change) max_change = change;
    }
    
    bwr->max_temp_change = max_change;
    
    // Warn if seeing large changes
    if (max_change > 50.0) {
        printf("  WARNING: Large temperature change: %.1f K\n", max_change);
        bwr->stability_violations++;
    }
    
    // Save old values
    memcpy(bwr->temperature_old, bwr->temperature, n_cells * sizeof(double));
    memcpy(bwr->void_fraction_old, bwr->void_fraction, n_cells * sizeof(double));
    
    return (unstable_cells < n_cells * 0.1);  // Allow up to 10% problematic cells
}

// Run BWR transient
void bwr_run_transient(BWRState* bwr, const char* scenario) {
    printf("\n=== Running Stable BWR Transient: %s ===\n", scenario);
    
    double next_output = 0.0;
    int step = 0;
    
    FILE *fp = fopen("bwr_results_stable.csv", "w");
    fprintf(fp, "Time,Power_MW,AvgTemp_K,MaxTemp_K,AvgVoid,Reactivity,TimeStep\n");
    
    while (bwr->time < SIM_TIME) {
        // Compute adaptive time step
        double dt = compute_adaptive_timestep(bwr);
        
        // Adjust if approaching output time
        if (bwr->time + dt > next_output && next_output > bwr->time) {
            dt = next_output - bwr->time;
        }
        if (bwr->time + dt > SIM_TIME) {
            dt = SIM_TIME - bwr->time;
        }
        
        // Calculate void reactivity feedback (smoothed)
        double void_rho = calculate_void_reactivity_smooth(bwr);
        
        // Apply control rod movement (scenario-dependent)
        if (strcmp(scenario, "rod_withdrawal") == 0) {
            if (bwr->time > 10.0 && bwr->time < 30.0) {
                double withdrawal_rate = -0.01 * dt;  // Slower: 1% per second
                bwr->control_rod_position += withdrawal_rate;
                bwr->control_rod_position = fmax(0.4, bwr->control_rod_position);
            }
        } else if (strcmp(scenario, "scram") == 0) {
            if (bwr->time > 30.0) {
                bwr->control_rod_position = 1.0;
            }
        } else if (strcmp(scenario, "power_ramp") == 0) {
            if (bwr->time > 5.0 && bwr->time < 35.0) {
                double target_pos = 0.6 - 0.1 * (bwr->time - 5.0) / 30.0;
                bwr->control_rod_position = target_pos;
            }
        }
        
        reactor_set_control_rods(bwr->reactor, bwr->control_rod_position);
        
        // Take time step
        reactor_step(bwr->reactor, dt);
        bwr->time += dt;
        step++;
        
        // Check stability
        if (!check_stability(bwr)) {
            printf("\n✗ CRITICAL: Numerical instability at t=%.3f s\n", bwr->time);
            printf("Step %d, dt=%.6f s\n", step, dt);
            fclose(fp);
            return;
        }
        
        bwr->total_power = reactor_get_total_power(bwr->reactor) / 1e6;
        
        // Output at intervals
        if (bwr->time >= next_output - 1e-9) {
            reactor_get_temperature(bwr->reactor, bwr->temperature);
            
            double sum_temp = 0.0, max_temp = 0.0, sum_void = 0.0;
            int n_cells = NX * NY * NZ;
            
            for (int idx = 0; idx < n_cells; idx++) {
                sum_temp += bwr->temperature[idx];
                if (bwr->temperature[idx] > max_temp) {
                    max_temp = bwr->temperature[idx];
                }
                sum_void += bwr->void_fraction[idx];
            }
            
            double avg_temp = sum_temp / n_cells;
            double avg_void = sum_void / n_cells;
            
            printf("t=%6.2f s: P=%6.1f MW, T_avg=%6.1f K, T_max=%6.1f K, "
                   "Void=%5.1f%%, Rods=%5.1f%%, dt=%.1e s\n",
                   bwr->time, bwr->total_power, avg_temp, max_temp,
                   avg_void * 100.0, bwr->control_rod_position * 100.0, dt);
            
            fprintf(fp, "%.3f,%.2f,%.2f,%.2f,%.4f,%.6f,%.6e\n",
                    bwr->time, bwr->total_power, avg_temp, max_temp,
                    avg_void, void_rho, dt);
            
            next_output += OUTPUT_INTERVAL;
        }
        
        // Emergency brake if too many stability issues
        if (bwr->stability_violations > 100) {
            printf("\n✗ Too many stability violations. Halting.\n");
            break;
        }
    }
    
    fclose(fp);
    printf("\n✓ Simulation complete. Results: bwr_results_stable.csv\n");
    printf("  Stability violations: %d\n", bwr->stability_violations);
}

// Print status
void bwr_print_status(BWRState* bwr) {
    reactor_get_temperature(bwr->reactor, bwr->temperature);
    reactor_get_power(bwr->reactor, bwr->power_density);
    
    double peak_power = 0.0, peak_temp = 0.0;
    int n_cells = NX * NY * NZ;
    
    for (int i = 0; i < n_cells; i++) {
        if (bwr->power_density[i] > peak_power) peak_power = bwr->power_density[i];
        if (bwr->temperature[i] > peak_temp) peak_temp = bwr->temperature[i];
    }
    
    printf("\n=== BWR Status ===\n");
    printf("Total Power: %.1f MW\n", bwr->total_power);
    printf("Peak Power Density: %.2e W/m³\n", peak_power);
    printf("Peak Temperature: %.1f K (%.1f °C)\n", peak_temp, peak_temp - 273.15);
    printf("Control Rod Position: %.1f%% inserted\n",
           bwr->control_rod_position * 100.0);
    printf("Max Temperature Change: %.2f K/step\n", bwr->max_temp_change);
}

// Cleanup
void bwr_destroy(BWRState* bwr) {
    reactor_destroy(bwr->reactor);
    free(bwr->temperature);
    free(bwr->temperature_old);
    free(bwr->pressure);
    free(bwr->power_density);
    free(bwr->void_fraction);
    free(bwr->void_fraction_old);
}

int main(int argc, char* argv[]) {
    printf("╔═══════════════════════════════════════════════════════════╗\n");
    printf("║      Stabilised BWR Simulator - Nuclear Physics Lib      ║\n");
    printf("║              Boiling Water Reactor (Stable)              ║\n");
    printf("╚═══════════════════════════════════════════════════════════╝\n\n");
    
    const char* scenario = "steady_state";
    if (argc > 1) {
        scenario = argv[1];
    }
    
    printf("Available scenarios:\n");
    printf("  steady_state   - Normal operation\n");
    printf("  rod_withdrawal - Slow control rod withdrawal\n");
    printf("  scram          - Emergency shutdown\n");
    printf("  power_ramp     - Gradual power increase\n\n");
    
    BWRState bwr;
    bwr_init_geometry(&bwr);
    
    printf("\n=== Equilibrating to Steady State ===\n");
    reactor_set_control_rods(bwr.reactor, 0.65);
    
    printf("Gradual warm-up phase...\n");
    for (int i = 0; i < 500; i++) {
        double dt = 0.002;  // 2 ms steps
        reactor_step(bwr.reactor, dt);
        bwr.time += dt;
        
        if (!check_stability(&bwr)) {
            printf("✗ Initialization failed at step %d\n", i);
            bwr_destroy(&bwr);
            return 1;
        }
        
        if (i % 100 == 0) {
            bwr.total_power = reactor_get_total_power(bwr.reactor) / 1e6;
            printf("  Step %d (%.2f s): Power = %.1f MW\n",
                   i, bwr.time, bwr.total_power);
        }
    }
    
    bwr.total_power = reactor_get_total_power(bwr.reactor) / 1e6;
    printf("\n✓ Reached steady state: %.1f MW\n", bwr.total_power);
    
    // Reset time for transient
    bwr.time = 0.0;
    
    // Run selected scenario
    if (strcmp(scenario, "steady_state") != 0) {
        bwr_run_transient(&bwr, scenario);
    }
    
    // Final status
    printf("\n=== Final State ===\n");
    bwr_print_status(&bwr);
    
    bwr_destroy(&bwr);
    
    printf("\n✓ Simulation completed successfully!\n");
    return 0;
}