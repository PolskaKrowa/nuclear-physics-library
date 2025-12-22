// examples/bwr_simulator.c
//
// Boiling Water Reactor (BWR) Simulator
//
// Simulates a simplified BWR with:
// - Core with fuel assemblies
// - Coolant flow and boiling
// - Steam separation
// - Control rod actuation
// - Void feedback reactivity
//
// Compile: gcc bwr_simulator.c -o bwr_sim -lnuclear_physics_c -lm
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nuclear_physics.h"

// BWR parameters
#define CORE_HEIGHT 3.7     // Core active height [m]
#define CORE_DIAMETER 4.5   // Core diameter [m]
#define FUEL_ASSEMBLIES 764 // Number of fuel assemblies
#define THERMAL_POWER 3293  // Thermal power [MW]

// Grid resolution
#define NX 30  // Radial
#define NY 30  // Radial
#define NZ 40  // Axial

// Physical constants
#define T_INLET 550.0   // Inlet temperature [K] (277°C)
#define T_SAT 560.0     // Saturation temperature [K] at ~7 MPa
#define P_OPERATING 7.0e6  // Operating pressure [Pa] (7 MPa)

// Simulation parameters
#define SIM_TIME 100.0     // Total simulation time [s]
#define OUTPUT_INTERVAL 1.0 // Output interval [s]

typedef struct {
    ReactorHandle reactor;
    int nx, ny, nz;
    double dx, dy, dz;
    double time;
    double total_power;
    double control_rod_position;
    
    // Arrays for monitoring
    double* temperature;
    double* pressure;
    double* power_density;
    double* void_fraction;
} BWRState;

// Initialize BWR geometry and materials
void bwr_init_geometry(BWRState* bwr) {
    int i, j, k;
    double r, x, y;
    double core_radius = CORE_DIAMETER / 2.0;
    
    printf("Initializing BWR geometry...\n");
    printf("  Core: %.2f m diameter x %.2f m height\n", CORE_DIAMETER, CORE_HEIGHT);
    printf("  Grid: %dx%dx%d (radial x radial x axial)\n", NX, NY, NZ);
    printf("  Thermal power: %.0f MW\n", THERMAL_POWER);
    
    // Create reactor
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
    bwr->control_rod_position = 0.0; // Fully withdrawn
    
    // Allocate monitoring arrays
    int n_cells = NX * NY * NZ;
    bwr->temperature = (double*)malloc(n_cells * sizeof(double));
    bwr->pressure = (double*)malloc(n_cells * sizeof(double));
    bwr->power_density = (double*)malloc(n_cells * sizeof(double));
    bwr->void_fraction = (double*)malloc(n_cells * sizeof(double));
    
    // Define fuel region (cylindrical core)
    double power_density = (THERMAL_POWER * 1e6) / 
                          (M_PI * core_radius * core_radius * CORE_HEIGHT);
    
    // Reduce by factor to account for actual enrichment and geometry
    power_density *= 0.1;  // Start lower for stability
    
    for (k = 1; k <= NZ; k++) {
        for (j = 1; j <= NY; j++) {
            for (i = 1; i <= NX; i++) {
                // Calculate radial position
                x = (i - NX/2.0) * bwr->dx;
                y = (j - NY/2.0) * bwr->dy;
                r = sqrt(x*x + y*y);
                
                // Core region (cylindrical)
                if (r < core_radius) {
                    // Fuel with radial and axial power peaking
                    double z_norm = (double)k / NZ;
                    double r_norm = r / core_radius;
                    
                    // Cosine axial profile, Bessel radial profile
                    double axial_factor = 1.0 + 0.3 * cos(M_PI * (z_norm - 0.5));
                    double radial_factor = 1.0 - 0.3 * r_norm * r_norm;
                    double local_power = power_density * axial_factor * radial_factor;
                    
                    reactor_set_fuel_region(bwr->reactor, i, i, j, j, k, k,
                                          local_power, 0.035); // 3.5% enriched
                } else {
                    // Reflector region (water)
                    reactor_set_coolant_region(bwr->reactor, i, i, j, j, k, k,
                                              T_INLET, 0.0);
                }
            }
        }
    }
    
    // Set coolant channels (interstitial spaces)
    reactor_set_coolant_region(bwr->reactor, 1, NX, 1, NY, 1, 5,
                               T_INLET, 15000.0); // Inlet plenum
    
    printf("Geometry initialized successfully.\n");
}

// Calculate void fraction from temperature
double calculate_void_fraction(double T) {
    if (T < T_SAT) {
        return 0.0; // Subcooled liquid
    } else if (T < T_SAT + 50.0) {
        // Linear transition to boiling
        return (T - T_SAT) / 50.0;
    } else {
        return 0.8; // Saturated steam-water mixture
    }
}

// Calculate void reactivity feedback
double calculate_void_reactivity(BWRState* bwr) {
    int idx, i, j, k;
    double total_void = 0.0;
    double core_radius = CORE_DIAMETER / 2.0;
    
    reactor_get_temperature(bwr->reactor, bwr->temperature);
    
    idx = 0;
    for (k = 0; k < NZ; k++) {
        for (j = 0; j < NY; j++) {
            for (i = 0; i < NX; i++) {
                double x = (i - NX/2.0) * bwr->dx;
                double y = (j - NY/2.0) * bwr->dy;
                double r = sqrt(x*x + y*y);
                
                if (r < core_radius) {
                    bwr->void_fraction[idx] = calculate_void_fraction(bwr->temperature[idx]);
                    total_void += bwr->void_fraction[idx];
                }
                idx++;
            }
        }
    }
    
    double avg_void = total_void / (NX * NY * NZ);
    
    // Void coefficient: typically -10 to -15 pcm/%void for BWR
    // Negative feedback: more void -> less moderation -> less reactivity
    double void_reactivity = -0.0015 * avg_void * 100.0; // Convert to absolute
    
    return void_reactivity;
}

// Run BWR transient simulation
void bwr_run_transient(BWRState* bwr, const char* scenario) {
    printf("\n=== Running BWR Transient: %s ===\n", scenario);
    
    double dt;
    double next_output = 0.0;
    int step = 0;
    
    FILE *fp = fopen("bwr_results.csv", "w");
    fprintf(fp, "Time,Power_MW,AvgTemp_K,MaxTemp_K,AvgVoid,Reactivity\n");
    
    while (bwr->time < SIM_TIME) {
        // Get adaptive time step
        dt = reactor_get_max_dt(bwr->reactor);
        dt = fmin(dt, 0.005);  // Limit to 5 ms for neutronics stability
        
        // Calculate void reactivity feedback
        double void_rho = calculate_void_reactivity(bwr);
        
        // Apply control rod movement (example scenarios)
        if (strcmp(scenario, "rod_withdrawal") == 0) {
            if (bwr->time > 10.0 && bwr->time < 30.0) {
                // Slow withdrawal over 20 seconds
                double withdrawal_rate = -0.02;  // 2% per second
                bwr->control_rod_position += withdrawal_rate * dt;
                bwr->control_rod_position = fmax(0.4, bwr->control_rod_position);
            }
        } else if (strcmp(scenario, "scram") == 0) {
            if (bwr->time > 30.0) {
                // Emergency shutdown (full insertion)
                bwr->control_rod_position = 1.0;
            }
        } else if (strcmp(scenario, "power_ramp") == 0) {
            if (bwr->time > 5.0 && bwr->time < 35.0) {
                // Slow power increase
                double target_pos = 0.6 - 0.15 * (bwr->time - 5.0) / 30.0;
                bwr->control_rod_position = target_pos;
            }
        }
        
        reactor_set_control_rods(bwr->reactor, bwr->control_rod_position);
        
        // Perform time step
        reactor_step(bwr->reactor, dt);
        bwr->time += dt;
        step++;
        
        // Safety check
        bwr->total_power = reactor_get_total_power(bwr->reactor) / 1e6;
        if (bwr->total_power > 1e6 || bwr->total_power != bwr->total_power) {
            printf("\n✗ ERROR: Numerical instability detected at t=%.2f s\n", bwr->time);
            printf("Power = %.2e MW (unrealistic)\n", bwr->total_power);
            fclose(fp);
            return;
        }
        
        // Output at intervals
        if (bwr->time >= next_output) {
            // Get current state
            reactor_get_temperature(bwr->reactor, bwr->temperature);
            
            // Calculate statistics
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
            
            printf("t=%.2f s: Power=%.1f MW, T_avg=%.1f K, T_max=%.1f K, "
                   "Void=%.1f%%, Rods=%.1f%%\n",
                   bwr->time, bwr->total_power, avg_temp, max_temp,
                   avg_void * 100.0, bwr->control_rod_position * 100.0);
            
            fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.4f,%.6f\n",
                    bwr->time, bwr->total_power, avg_temp, max_temp,
                    avg_void, void_rho);
            
            next_output += OUTPUT_INTERVAL;
        }
    }
    
    fclose(fp);
    printf("\nSimulation complete. Results saved to bwr_results.csv\n");
}

// Print BWR status
void bwr_print_status(BWRState* bwr) {
    reactor_get_temperature(bwr->reactor, bwr->temperature);
    reactor_get_pressure(bwr->reactor, bwr->pressure);
    reactor_get_power(bwr->reactor, bwr->power_density);
    
    // Find peak values
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
    printf("Control Rod Position: %.1f%% inserted\n", bwr->control_rod_position * 100.0);
}

// Clean up BWR simulation
void bwr_destroy(BWRState* bwr) {
    reactor_destroy(bwr->reactor);
    free(bwr->temperature);
    free(bwr->pressure);
    free(bwr->power_density);
    free(bwr->void_fraction);
}

int main(int argc, char* argv[]) {
    printf("╔═══════════════════════════════════════════════════════════╗\n");
    printf("║         BWR Simulator - Nuclear Physics Library          ║\n");
    printf("║                 Boiling Water Reactor                     ║\n");
    printf("╚═══════════════════════════════════════════════════════════╝\n\n");
    
    // Parse command line
    const char* scenario = "steady_state";
    if (argc > 1) {
        scenario = argv[1];
    }
    
    printf("Available scenarios:\n");
    printf("  steady_state   - Normal operation\n");
    printf("  rod_withdrawal - Control rod withdrawal transient\n");
    printf("  scram          - Emergency shutdown\n");
    printf("  power_ramp     - Slow power increase\n\n");
    
    // Initialize BWR
    BWRState bwr;
    bwr_init_geometry(&bwr);
    
    // Run steady state first
    printf("\n=== Initializing to Steady State ===\n");
    
    // Start subcritical with rods partially inserted
    reactor_set_control_rods(bwr.reactor, 0.85);  // 85% inserted = subcritical
    
    printf("Phase 1: Subcritical hold (rods 85%% inserted)...\n");
    for (int i = 0; i < 50; i++) {
        double dt = 0.001;  // 1 ms steps for stability
        reactor_step(bwr.reactor, dt);
        bwr.time += dt;
    }
    
    // Gradually approach critical
    printf("Phase 2: Approaching critical...\n");
    for (int i = 0; i < 200; i++) {
        double rod_pos = 0.85 - (double)i / 200.0 * 0.25;  // Withdraw to 60%
        reactor_set_control_rods(bwr.reactor, rod_pos);
        
        double dt = 0.001;
        reactor_step(bwr.reactor, dt);
        bwr.time += dt;
        
        if (i % 50 == 0) {
            bwr.total_power = reactor_get_total_power(bwr.reactor) / 1e6;
            printf("  Step %d: Power = %.1f MW, Rods = %.1f%%\n", 
                   i, bwr.total_power, rod_pos * 100.0);
        }
    }
    
    // Hold at critical for equilibration
    printf("Phase 3: Critical equilibration...\n");
    reactor_set_control_rods(bwr.reactor, 0.60);  // 60% for critical
    
    for (int i = 0; i < 500; i++) {
        double dt = 0.01;
        reactor_step(bwr.reactor, dt);
        bwr.time += dt;
        
        if (i % 100 == 0) {
            bwr.total_power = reactor_get_total_power(bwr.reactor) / 1e6;
            printf("  Power = %.1f MW\n", bwr.total_power);
        }
    }
    
    bwr.total_power = reactor_get_total_power(bwr.reactor) / 1e6;
    
    if (bwr.total_power > 1e10 || bwr.total_power != bwr.total_power) {
        printf("\n✗ ERROR: Reactor went supercritical! Power = %.2e MW\n", bwr.total_power);
        printf("This indicates numerical instability.\n");
        printf("Try reducing time step or adjusting reactivity.\n");
        bwr_destroy(&bwr);
        return 1;
    }
    
    printf("\n✓ Reached steady state: %.1f MW\n", bwr.total_power);
    bwr.control_rod_position = 0.60;
    
    // Reset time for transient
    bwr.time = 0.0;
    
    // Run selected transient
    if (strcmp(scenario, "steady_state") != 0) {
        bwr_run_transient(&bwr, scenario);
    } else {
        bwr_print_status(&bwr);
    }
    
    // Final status
    printf("\n=== Final State ===\n");
    bwr_print_status(&bwr);
    
    // Cleanup
    bwr_destroy(&bwr);
    
    printf("\nSimulation complete!\n");
    return 0;
}