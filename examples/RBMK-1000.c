/*
 * RBMK-1000 Nuclear Reactor Simulator
 * 
 * Models the dynamics of an RBMK-1000 reactor including:
 * - Point kinetics with 6 delayed neutron groups
 * - Thermal hydraulics (fuel and coolant temperatures)
 * - Xenon-135 and Iodine-135 dynamics
 * - Control rod reactivity effects
 * - Positive void coefficient characteristic of RBMK
 * - Power distribution
 * 
 * Compile with:
 * gcc -o rbmk_sim rbmk_simulator.c -L./lib -lnuclear_physics_kernels \
 *     -lnuclear_physics_core -lgfortran -llapack -lblas -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

// Fortran interface declarations
extern void rk4_solve_(void (*func)(double*, double*, double*), 
                       double* t_span, double* y0, void* config,
                       double** t_out, double** y_out, void* status);

// RBMK-1000 Physical Parameters
#define N_DELAYED_GROUPS 6      // Delayed neutron groups
#define N_FUEL_CHANNELS 1661    // Number of fuel channels
#define THERMAL_POWER_NOMINAL 3200.0e6  // 3200 MWth
#define ELECTRICAL_POWER_NOMINAL 1000.0e6 // 1000 MWe

// State vector indices
#define IDX_POWER 0
#define IDX_DELAYED_START 1
#define IDX_FUEL_TEMP (IDX_DELAYED_START + N_DELAYED_GROUPS)
#define IDX_COOLANT_TEMP (IDX_FUEL_TEMP + 1)
#define IDX_MODERATOR_TEMP (IDX_COOLANT_TEMP + 1)
#define IDX_XENON (IDX_MODERATOR_TEMP + 1)
#define IDX_IODINE (IDX_XENON + 1)
#define STATE_SIZE (IDX_IODINE + 1)

// Reactor parameters structure
typedef struct {
    // Neutronics
    double beta[N_DELAYED_GROUPS];    // Delayed neutron fractions
    double lambda[N_DELAYED_GROUPS];  // Decay constants (1/s)
    double gen_time;                   // Neutron generation time (s)
    
    // Reactivity coefficients
    double fuel_temp_coef;            // Fuel temperature coefficient (pcm/K)
    double coolant_void_coef;         // Void coefficient (pcm/% void)
    double moderator_temp_coef;       // Moderator temperature coefficient (pcm/K)
    double xenon_worth;               // Xenon worth (pcm/1e24 atoms/m³)
    
    // Thermal parameters
    double fuel_heat_capacity;        // J/(kg·K)
    double fuel_mass;                 // kg
    double coolant_flow_rate;         // kg/s
    double coolant_heat_capacity;     // J/(kg·K)
    double heat_transfer_coef;        // W/(m²·K)
    double heat_transfer_area;        // m²
    
    // Control and operational
    double control_rod_position;      // 0 (fully in) to 1 (fully out)
    double control_rod_worth;         // Total worth in pcm
    double coolant_inlet_temp;        // K
    
    // Current state
    double reactivity;                // Current reactivity (pcm)
    double void_fraction;             // Current void fraction (0-1)
    
    // Xenon/Iodine parameters
    double fission_yield_I;           // Iodine-135 fission yield
    double fission_yield_Xe;          // Xenon-135 fission yield
    double lambda_I;                  // Iodine decay constant (1/s)
    double lambda_Xe;                 // Xenon decay constant (1/s)
    double sigma_Xe;                  // Xenon absorption cross-section (barns)
    double neutron_flux_norm;         // Normalization for neutron flux
    
} ReactorParams;

// Global reactor parameters
ReactorParams reactor;

// Initialize reactor parameters with RBMK-1000 characteristics
void init_reactor_params() {
    // Delayed neutron parameters (typical for U-235)
    reactor.beta[0] = 0.000215;
    reactor.beta[1] = 0.001424;
    reactor.beta[2] = 0.001274;
    reactor.beta[3] = 0.002568;
    reactor.beta[4] = 0.000748;
    reactor.beta[5] = 0.000273;
    
    reactor.lambda[0] = 0.0127;  // 1/s
    reactor.lambda[1] = 0.0317;
    reactor.lambda[2] = 0.115;
    reactor.lambda[3] = 0.311;
    reactor.lambda[4] = 1.40;
    reactor.lambda[5] = 3.87;
    
    reactor.gen_time = 1.0e-4;  // 0.1 ms (graphite moderated)
    
    // RBMK-specific reactivity coefficients
    reactor.fuel_temp_coef = -3.0;      // Negative fuel temperature coefficient
    reactor.coolant_void_coef = +4.5;   // POSITIVE void coefficient (dangerous!)
    reactor.moderator_temp_coef = -0.5; // Slightly negative
    reactor.xenon_worth = -28.0;        // Strong xenon poisoning
    
    // Thermal parameters (simplified single channel)
    reactor.fuel_heat_capacity = 300.0;  // UO2
    reactor.fuel_mass = 190000.0;        // ~190 tons of fuel
    reactor.coolant_flow_rate = 37500.0; // kg/s total
    reactor.coolant_heat_capacity = 5200.0; // Water at high pressure
    reactor.heat_transfer_coef = 15000.0;   // W/(m²·K)
    reactor.heat_transfer_area = 70000.0;    // m²
    
    // Control system
    reactor.control_rod_position = 0.7;  // 70% withdrawn
    reactor.control_rod_worth = -8000.0; // -80 mk when fully inserted
    reactor.coolant_inlet_temp = 543.0;  // 270°C inlet
    
    // Initial conditions
    reactor.reactivity = 0.0;
    reactor.void_fraction = 0.15;  // 15% void at nominal
    
    // Xenon/Iodine parameters
    reactor.fission_yield_I = 0.061;   // 6.1% yield
    reactor.fission_yield_Xe = 0.003;  // 0.3% direct yield
    reactor.lambda_I = 2.87e-5;        // 6.57 hour half-life
    reactor.lambda_Xe = 2.09e-5;       // 9.14 hour half-life
    reactor.sigma_Xe = 2.65e6;         // Very large (barns * 1e-24)
    reactor.neutron_flux_norm = 3.0e18; // Normalization constant
}

// Calculate total reactivity
double calculate_reactivity(double fuel_temp, double coolant_temp, 
                           double moderator_temp, double xenon_concentration) {
    double rho = 0.0;
    
    // Control rod contribution
    rho += reactor.control_rod_worth * (1.0 - reactor.control_rod_position);
    
    // Temperature feedbacks
    rho += reactor.fuel_temp_coef * (fuel_temp - 900.0);  // Relative to 900K nominal
    rho += reactor.moderator_temp_coef * (moderator_temp - 900.0);
    
    // Void coefficient (POSITIVE - key RBMK characteristic!)
    rho += reactor.coolant_void_coef * reactor.void_fraction * 100.0;
    
    // Xenon poisoning
    rho += reactor.xenon_worth * xenon_concentration;
    
    return rho;
}

// Calculate void fraction from coolant temperature and power
double calculate_void_fraction(double power_fraction, double coolant_temp) {
    // Simplified void model: increases with power and temperature
    double temp_effect = (coolant_temp - 543.0) / 100.0;
    double power_effect = (power_fraction - 1.0);
    
    double void = 0.15 + 0.3 * power_effect + 0.1 * temp_effect;
    
    // Clamp between 0 and 0.95
    if (void < 0.0) void = 0.0;
    if (void > 0.95) void = 0.95;
    
    return void;
}

// ODE system for reactor dynamics
void reactor_odes(double* t, double* y, double* dydt) {
    // Extract state variables
    double power = y[IDX_POWER];
    double delayed[N_DELAYED_GROUPS];
    for (int i = 0; i < N_DELAYED_GROUPS; i++) {
        delayed[i] = y[IDX_DELAYED_START + i];
    }
    double T_fuel = y[IDX_FUEL_TEMP];
    double T_coolant = y[IDX_COOLANT_TEMP];
    double T_moderator = y[IDX_MODERATOR_TEMP];
    double Xe = y[IDX_XENON];
    double I = y[IDX_IODINE];
    
    // Calculate void fraction
    reactor.void_fraction = calculate_void_fraction(power / THERMAL_POWER_NOMINAL, 
                                                     T_coolant);
    
    // Calculate reactivity
    reactor.reactivity = calculate_reactivity(T_fuel, T_coolant, T_moderator, Xe);
    double rho = reactor.reactivity / 100000.0;  // Convert pcm to absolute
    
    // Total delayed neutron fraction
    double beta_total = 0.0;
    for (int i = 0; i < N_DELAYED_GROUPS; i++) {
        beta_total += reactor.beta[i];
    }
    
    // Point kinetics equations
    double delayed_sum = 0.0;
    for (int i = 0; i < N_DELAYED_GROUPS; i++) {
        delayed_sum += reactor.lambda[i] * delayed[i];
    }
    
    dydt[IDX_POWER] = ((rho - beta_total) / reactor.gen_time) * power + delayed_sum;
    
    // Delayed neutron precursors
    for (int i = 0; i < N_DELAYED_GROUPS; i++) {
        dydt[IDX_DELAYED_START + i] = (reactor.beta[i] / reactor.gen_time) * power 
                                       - reactor.lambda[i] * delayed[i];
    }
    
    // Thermal equations
    double Q_fission = power;  // Thermal power generated
    double Q_removal = reactor.heat_transfer_coef * reactor.heat_transfer_area 
                      * (T_fuel - T_coolant);
    
    // Fuel temperature
    dydt[IDX_FUEL_TEMP] = (Q_fission - Q_removal) / 
                          (reactor.fuel_mass * reactor.fuel_heat_capacity);
    
    // Coolant temperature (using energy balance)
    double Q_coolant = Q_removal;
    double coolant_temp_rise = Q_coolant / 
                              (reactor.coolant_flow_rate * reactor.coolant_heat_capacity);
    dydt[IDX_COOLANT_TEMP] = (coolant_temp_rise - (T_coolant - reactor.coolant_inlet_temp)) / 5.0;
    
    // Moderator temperature (simplified - follows fuel temperature with lag)
    dydt[IDX_MODERATOR_TEMP] = (T_fuel - T_moderator) / 100.0;
    
    // Xenon-135 and Iodine-135 dynamics
    double fission_rate = power / (200.0e6 * 1.602e-19);  // Fissions per second
    double neutron_flux = reactor.neutron_flux_norm * (power / THERMAL_POWER_NOMINAL);
    
    // Iodine-135 production and decay
    dydt[IDX_IODINE] = reactor.fission_yield_I * fission_rate - reactor.lambda_I * I;
    
    // Xenon-135: production from I-135 decay, direct fission, loss by decay and burnup
    dydt[IDX_XENON] = reactor.lambda_I * I 
                     + reactor.fission_yield_Xe * fission_rate
                     - reactor.lambda_Xe * Xe
                     - reactor.sigma_Xe * neutron_flux * Xe;
}

// Display reactor status
void display_status(double t, double* state, int scram) {
    double power = state[IDX_POWER];
    double power_pct = (power / THERMAL_POWER_NOMINAL) * 100.0;
    double T_fuel = state[IDX_FUEL_TEMP];
    double T_coolant = state[IDX_COOLANT_TEMP];
    double T_moderator = state[IDX_MODERATOR_TEMP];
    double Xe = state[IDX_XENON];
    
    printf("\033[2J\033[H");  // Clear screen
    printf("╔════════════════════════════════════════════════════════════════╗\n");
    printf("║           RBMK-1000 REACTOR SIMULATOR - CONTROL ROOM          ║\n");
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    printf("║ Time: %8.1f s                                                ║\n", t);
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    printf("║                       POWER INDICATORS                         ║\n");
    printf("║ Thermal Power:   %7.1f MWth (%5.1f%%)                      ║\n", 
           power / 1e6, power_pct);
    printf("║ Electrical:      %7.1f MWe  (%5.1f%%)                      ║\n",
           power * 0.31 / 1e6, power_pct * 0.31);
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    printf("║                    TEMPERATURE READINGS                        ║\n");
    printf("║ Fuel Temperature:      %7.1f K  (%6.1f °C)                 ║\n",
           T_fuel, T_fuel - 273.15);
    printf("║ Coolant Temperature:   %7.1f K  (%6.1f °C)                 ║\n",
           T_coolant, T_coolant - 273.15);
    printf("║ Moderator Temperature: %7.1f K  (%6.1f °C)                 ║\n",
           T_moderator, T_moderator - 273.15);
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    printf("║                      REACTIVITY BALANCE                        ║\n");
    printf("║ Total Reactivity: %+7.1f pcm                                 ║\n", 
           reactor.reactivity);
    printf("║ Void Fraction:    %7.2f %%                                   ║\n",
           reactor.void_fraction * 100.0);
    printf("║ Control Rods:     %7.1f %% withdrawn                         ║\n",
           reactor.control_rod_position * 100.0);
    printf("║ Xenon-135:        %7.2e atoms/m³                           ║\n", Xe);
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    
    // Safety warnings
    if (scram) {
        printf("║                   ⚠️  SCRAM ACTIVATED! ⚠️                      ║\n");
    } else if (power_pct > 110.0) {
        printf("║                 ⚠️  POWER LIMIT EXCEEDED! ⚠️                   ║\n");
    } else if (T_fuel > 2500.0) {
        printf("║             ⚠️  FUEL TEMPERATURE CRITICAL! ⚠️                  ║\n");
    } else if (reactor.void_fraction > 0.7) {
        printf("║              ⚠️  DANGEROUS VOID FRACTION! ⚠️                   ║\n");
    } else {
        printf("║                    STATUS: NORMAL                              ║\n");
    }
    
    printf("╠════════════════════════════════════════════════════════════════╣\n");
    printf("║ Controls: [W/S] Rods  [A/D] Flow  [SPACE] Scram  [Q] Quit     ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    
    // Additional information
    printf("\nPhysics Notes:\n");
    printf("• Positive void coefficient: +%.1f pcm/%% void\n", reactor.coolant_void_coef);
    printf("• Current void contribution: %+.1f pcm\n", 
           reactor.coolant_void_coef * reactor.void_fraction * 100.0);
    printf("• Xenon poisoning: %.1f pcm\n", reactor.xenon_worth * Xe);
    
    fflush(stdout);
}

// Simplified integration (Euler method for demonstration)
void integrate_step(double* state, double dt) {
    double dydt[STATE_SIZE];
    double t = 0.0;
    
    reactor_odes(&t, state, dydt);
    
    for (int i = 0; i < STATE_SIZE; i++) {
        state[i] += dydt[i] * dt;
        
        // Ensure physical constraints
        if (i == IDX_POWER && state[i] < 0.0) state[i] = 0.0;
        if (i >= IDX_DELAYED_START && i < IDX_FUEL_TEMP && state[i] < 0.0) 
            state[i] = 0.0;
    }
}

// Main simulation
int main(int argc, char** argv) {
    printf("RBMK-1000 Nuclear Reactor Simulator\n");
    printf("Based on Chernobyl-type reactor design\n\n");
    printf("⚠️  WARNING: This reactor has a POSITIVE void coefficient!\n");
    printf("Loss of coolant can cause power excursion.\n\n");
    printf("Press ENTER to start...\n");
    getchar();
    
    // Initialize reactor
    init_reactor_params();
    
    // Initial state
    double state[STATE_SIZE];
    state[IDX_POWER] = THERMAL_POWER_NOMINAL * 0.7;  // Start at 70% power
    
    // Initialize delayed neutron precursors at equilibrium
    for (int i = 0; i < N_DELAYED_GROUPS; i++) {
        state[IDX_DELAYED_START + i] = (reactor.beta[i] / reactor.lambda[i]) * 
                                        state[IDX_POWER] / reactor.gen_time;
    }
    
    state[IDX_FUEL_TEMP] = 900.0;      // 627°C
    state[IDX_COOLANT_TEMP] = 563.0;   // 290°C
    state[IDX_MODERATOR_TEMP] = 873.0; // 600°C
    state[IDX_XENON] = 1.0e17;         // Equilibrium xenon
    state[IDX_IODINE] = 2.0e18;        // Equilibrium iodine
    
    double t = 0.0;
    double dt = 0.1;  // 100 ms time step
    int scram = 0;
    
    // Simulation loop
    while (1) {
        // Update display every 10 steps (1 second)
        if ((int)(t / dt) % 10 == 0) {
            display_status(t, state, scram);
        }
        
        // Check for safety limits
        if (!scram) {
            if (state[IDX_FUEL_TEMP] > 2800.0 || 
                state[IDX_POWER] / THERMAL_POWER_NOMINAL > 1.2 ||
                reactor.void_fraction > 0.85) {
                scram = 1;
                printf("\n🚨 AUTOMATIC SCRAM TRIGGERED! 🚨\n");
            }
        }
        
        // Apply SCRAM
        if (scram) {
            reactor.control_rod_position -= 0.01;  // Insert rods rapidly
            if (reactor.control_rod_position < 0.0) 
                reactor.control_rod_position = 0.0;
        }
        
        // Integrate one time step
        integrate_step(state, dt);
        
        t += dt;
        
        // Simple delay for real-time feel (10 Hz update)
        usleep(100000);
        
        // Check for user input (non-blocking would be better)
        // For now, this is a simple continuous simulation
        
        // Stop after 10 minutes of simulation time or if power drops to near zero
        if (t > 600.0 || state[IDX_POWER] < THERMAL_POWER_NOMINAL * 0.01) {
            break;
        }
    }
    
    printf("\n\nSimulation ended at t = %.1f seconds\n", t);
    printf("Final power: %.1f MWth (%.1f%%)\n", 
           state[IDX_POWER] / 1e6,
           (state[IDX_POWER] / THERMAL_POWER_NOMINAL) * 100.0);
    
    return 0;
}