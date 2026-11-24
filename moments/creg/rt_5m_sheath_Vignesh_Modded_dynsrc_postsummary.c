#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <gkyl_const.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h> 
#endif

#include <rt_arg_parse.h>

// heresy (this is a private .h file I shouldn't really be using it but I want to access values that require app to not be an opaque pointer)
#include <gkyl_moment_priv.h>

struct profile_data {
    int n_points;
    double *data;
};

struct sheath_ctx
{
    // Mathematical constants (dimensionless).
    double pi;

    // Physical constants (using normalized code units).
    double gas_Gamma; // Adiabatic index.
    double ep0;       // Permittivity of free space.
    double mu0;       // Permeability of free space.

    // Electron species parameters.
    double me;  // Electron mass.
    double qe;  // Electron charge.
    double ne;  // Electron number density.
    double Te;  // Electron temperature.
    double vte; // Electron thermal velocity.
    double ue;  // Electron drift velocity.

    // Ion species parameters.
    double mi;  // Ion mass.
    double qi;  // Ion charge.
    double ni;  // Ion number density.
    double Ti;  // Ion temperature.
    double vti; // Ion thermal velocity.
    double ui;  // Ion drift velocity.

    // Derived physical quantities (using normalized code units).
    double wpe;     // Electron plasma frequency.
    double lambdaD; // Debye length.
    double ub;      // Bohm speed.

    // Collision parameters.
    bool has_collision;
    double nu_base_ei;

    // Simulation parameters.
    double Lx;
    int Nx;
    double cfl_frac;

    double t_end;
    int num_frames;
    int field_energy_calcs;   // Number of times to calculate field energy.
    int integrated_mom_calcs; // Number of times to calculate integrated moments.
    double dt_failure_tol;
    int num_failures_max;

    struct profile_data data_ne;
    struct profile_data data_ni;
    struct profile_data data_Te;
    struct profile_data data_Ti;
    struct profile_data data_ui;

    // Flux storage variables (for each time step)
    double out_e_flux;
    double out_i_flux;
    double in_e_flux;
    double in_i_flux;

    // Pre-sheath length
    double presheath_length;
    double source_domain_length;
    double transition_length;

    // Particle count variables
    double e_count;
    double i_count;

    // Runtime sheath profile structs
    struct profile_data elc_sheath;
    struct profile_data ion_sheath;
};

// Function to read data from file
void read_profile(const char *fname, struct profile_data *p) {
    FILE *fp = fopen(fname, "r");
    if (!fp) {
        fprintf(stderr, "Error opening file %s\n", fname);
        exit(1);
    }

    // First pass: count points
    int count = 0;
    double val;
    while (fscanf(fp, "%lf", &val) == 1) {
        count++;
    }
    rewind(fp);

    p->n_points = count;
    p->data = (double *)gkyl_malloc(count * sizeof(double));

    // Second pass: read data
    int i = 0;
    while (fscanf(fp, "%lf", &p->data[i]) == 1) {
        i++;
    }
    fclose(fp);
    printf("Read %d points from %s\n", count, fname);
}

// Function to interpolate value at x
// Assumes data covers range [0, Lx] uniformly for (int i = 0; i < NX; i++) {
double get_interp_val(struct profile_data *p, double x, double Lx) {
    if (p->n_points == 0) return 0.0;
    
    // Normalize x to [0, 1] then map to index
    double x_norm = x / Lx;
    if (x_norm < 0.0) x_norm = 0.0;
    if (x_norm > 1.0) x_norm = 1.0;

    double idx_f = x_norm * (p->n_points - 1);
    int idx = (int)idx_f;
    
    if (idx >= p->n_points - 1) return p->data[p->n_points - 1];

    double frac = idx_f - idx;
    return (1.0 - frac) * p->data[idx] + frac * p->data[idx + 1];
}

struct sheath_ctx
create_ctx(void)
{
     
    // Mathematical constants (dimensionless).
    double pi = M_PI;

    // Physical constants (using normalized code units).
    double gas_Gamma = 5.0 / 3.0; // Adiabatic index.
    double ep0 = 1.0;             // Permittivity of free space.
    double mu0 = 0.01;       // Permeability of free space.

    // Electron species parameters.
    double me = 1.0;            // Electron mass.
    double qe = -1.0;           // Electron charge.
    double ne = 1.0;            // Electron number density.
    double Te = 1.0;           // Electron temperature.
    double vte = sqrt(Te / me); // Electron thermal velocity.
    double ue = 0.0 * vte;      // Electron drift velocity.

    // Ion species parameters.
    double mi = 1836.2 * me;    // Ion mass.
    double qi = -qe;            // Ion charge.
    double ni = ne;            // Ion number density.
    double Ti = Te;            // Ion temperature.
    double vti = sqrt(Ti / mi); // Ion thermal velocity.
    double ui = 0.0 * vti;      // Ion drift velocity.

    // Derived physical quantities (using normalized code units).
    double wpe = sqrt(ne * qe * qe / (me * ep0)); // Electron plasma frequency.
    double lambdaD = vte / wpe;                   // Debye length.
    double ub = sqrt((Te + Ti) / mi);              // Bohm speed.

    // Collision parameters.
    bool has_collision = false;
    double mfp = 5.0 * lambdaD;
    double nu_base_ei = (sqrt(2.0) * vte) / mfp;

    // Simulation parameters.
    double Lx = lambdaD * 128.0; // Domain size (x-direction).
    int Nx = Lx*2;                  // Cell count (x-direction). Original was Lx * 2
    double cfl_frac = 0.95;            // CFL coefficient. Original was 0.95

    double t_end = 30000.0 / wpe;   // Final simulation time.
    int num_frames = 100;            // Number of output frames.
    int field_energy_calcs = INT_MAX;   // Number of times to calculate field energy.
    int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
    double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
    int num_failures_max = 20;      // Maximum allowable number of consecutive small time-steps.

    // Initializing pre-sheath length
    double presheath_length = 50.0;
    double source_domain_length = 50.0;
    double transition_length = 16.0;
    
    struct sheath_ctx ctx = {
        .pi = pi,

        .gas_Gamma = gas_Gamma,
        .ep0 = ep0,
        .mu0 = mu0, 

        .me = me,
        .qe = qe,
        .ne = ne,
        .Te = Te,
        .vte = vte,
        .ue = ue,

        .mi = mi,
        .qi = qi,
        .ni = ni,
        .Ti = Ti,
        .vti = vti,
        .ui = ui,

        .wpe = wpe,
        .lambdaD = lambdaD,
        .ub = ub,

        .has_collision = has_collision,
        .nu_base_ei = nu_base_ei,

        .Lx = Lx,
        .Nx = Nx,
        .cfl_frac = cfl_frac,

        .t_end = t_end,
        .num_frames = num_frames,
        .field_energy_calcs = field_energy_calcs,
        .integrated_mom_calcs = integrated_mom_calcs,
        .dt_failure_tol = dt_failure_tol,
        .num_failures_max = num_failures_max,

        .presheath_length = presheath_length,
        .source_domain_length = source_domain_length,
        .transition_length = transition_length,
       
    };

// Load Profile Data
    // Assumes files are in the local directory
    read_profile("ne_Rn.txt", &ctx.data_ne);
    read_profile("ni_Rn.txt", &ctx.data_ni);
    read_profile("Te_Rn.txt", &ctx.data_Te);
    read_profile("Ti_Rn.txt", &ctx.data_Ti);
    read_profile("ui_Rn.txt", &ctx.data_ui);

    // Initializing net flux values
    double out_e_flux = 0.0;
    double out_i_flux = 0.0;
    double in_e_flux = 0.0;
    double in_i_flux = 0.0;

    // Initializing particle count variables
    double e_count = 0.0;
    double i_count = 0.0;

    // Initializing runtime sheath profiles
    ctx.elc_sheath.n_points = Nx;
    ctx.ion_sheath.n_points = Nx;
    ctx.elc_sheath.data = (double *)gkyl_malloc(Nx * sizeof(double));
    ctx.ion_sheath.data = (double *)gkyl_malloc(Nx * sizeof(double));
    for (int i = 0; i < Nx; i++) {
        ctx.elc_sheath.data[i] = 0.0;
        ctx.ion_sheath.data[i] = 0.0;
    }
return ctx;
}
void evalElcInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
    struct sheath_ctx *app = ctx;
    double x = xn[0];

    // Interpolate values from loaded data
    // Assuming ne.txt is normalized to n0 (app->ne)
    double ne_norm = get_interp_val(&app->data_ne, x, app->Lx);
    // Assuming Te.txt is in Joules (SI) based on magnitude ~1e-17
    double Te_si = get_interp_val(&app->data_Te, x, app->Lx);

    double n_local = ne_norm * app->ne;
    
    // Safety check for vacuum regions if any
    if (n_local < 1.0e-5 * app->ne) n_local = 1.0e-5 * app->ne;

    double rhoe = n_local * app->me;
    double Ee = n_local * Te_si / (app->gas_Gamma - 1.0);

    // Set electron mass density.
    fout[0] = rhoe;
    // Set electron momentum density.
    fout[1] = 0.0;
    fout[2] = 0.0;
    fout[3] = 0.0;
    // Set electron total energy density.
    fout[4] = Ee;
}

void evalIonInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
    struct sheath_ctx *app = ctx;
    double x = xn[0];

    // Interpolate values from loaded data
    double ni_norm = get_interp_val(&app->data_ni, x, app->Lx);
    double Ti_si = get_interp_val(&app->data_Ti, x, app->Lx);
    double ui_norm = get_interp_val(&app->data_ui, x, app->Lx);

    double n_local = ni_norm * app->ni;
    
    // Scale velocity: assuming ui in file is normalized to Bohm speed (ub)
    // If ui in file is raw SI (m/s), remove (* app->ub)
    double u_local = ui_norm * app->ub;

    if (n_local < 1.0e-5 * app->ni) n_local = 1.0e-5 * app->ni;

    double rhoi = n_local * app->mi;
    double momi_x = rhoi * u_local;
    double Ei = n_local * Ti_si / (app->gas_Gamma - 1.0) + 0.5 * rhoi * u_local * u_local;

    // Set ion mass density.
    fout[0] = rhoi;
    // Set ion momentum density.
    fout[1] = momi_x;
    fout[2] = 0.0;
    fout[3] = 0.0;
    // Set ion total energy density.
    fout[4] = Ei;
}

void evalFieldInit(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
    struct sheath_ctx *app = ctx;
    double x = xn[0];
    // Set electric field.
    fout[0] = 0.0, fout[1] = 0.0, fout[2] = 0.0;
    // Set magnetic field.
    fout[3] = 0.0, fout[4] = 0.0, fout[5] = 0.0;
    // Set correction potentials.
    fout[6] = 0.0, fout[7] = 0.0;
}
// Sheath Electron
void evalElcUpperBC(const struct gkyl_wv_eqn *eqn, double t, int nc, const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
    struct sheath_ctx *app = ctx;

    double rhoe_BC = 1.0 * app->ne * app->me *  1e-19;
    double Ee_BC = 1.0 * app->ne * app->Te  * 1e-19 / (app->gas_Gamma - 1.0) ;

    // printf("Elc upper ghost cell 1 value (pre-write): %7.7f\n", ghost[1]);
    double ghost_mom = ghost[1];

    // Set Elc upper  boundary mass density.
    ghost[0] = rhoe_BC;
    // Set Elc upper  boundary momentum density.
    ghost[1] = 0.0;
    ghost[2] = 0.0;
    ghost[3] = 0.0;
    // Set Elc upper boundary total energy density.
    ghost[4] = Ee_BC;

    app->out_e_flux = skin[1];
}
// Source Electron
void evalElcLowerBC(const struct gkyl_wv_eqn *eqn, double t, int nc, const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
    struct sheath_ctx *app = ctx;

    double rhoe_BC = 1 * app->ne * app->me;
    double Ee_BC = 1 * app->ne * app->Te / (app->gas_Gamma - 1.0);

    // try printing out the values in the ghost cells before overwriting
    // Also, see if reported flux actually equals lost particles
    double ghost_mom = ghost[1];

    // Set Elc lower boundary mass density.
    /* ghost[0] = rhoe_BC;
    // Set Elc lower boundary momentum density.
    ghost[1] = 0.0;
    ghost[2] = 0.0;
    ghost[3] = 0.0;
    // Set Elc lower boundary total energy density.
    ghost[4] = Ee_BC; */

    // Copying the GKYL_SPECIES_COPY boundary condition requirements
    ghost[0] = skin[0];
    ghost[1] = skin[1];
    ghost[2] = skin[2];
    ghost[3] = skin[3];
    ghost[4] = skin[4];

    // Testing with only momentum at 0.0
    ghost[1] = 0.0;
    ghost[2] = 0.0;
    ghost[3] = 0.0;

    app->in_e_flux = skin[1];
    // printf("Elc in-flux w/ difference: %7.7f\n", skin[1]);
}
//Source Ion
void evalIonLowerBC(const struct gkyl_wv_eqn *eqn, double t, int nc, const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
    struct sheath_ctx *app = ctx;

    double rhoi_BC = 1.0 * app->ni * app->mi;
    double Ei_BC = 1.0 * app->ni * app->Ti / (app->gas_Gamma - 1.0);
    double momi_x = rhoi_BC * app->ui;

    double ghost_mom = ghost[1];

    // Set ion lower boundary mass density.
    /* ghost[0] = rhoi_BC;
    // Set ion lower boundary momentum density.
    ghost[1] = momi_x;
    ghost[2] = 0.0;
    ghost[3] = 0.0;
    // Set ion lower boundary total energy density.
    ghost[4] = Ei_BC; */

    // Copying the GKYL_SPECIES_COPY boundary condition requirements
    ghost[0] = skin[0];
    ghost[1] = skin[1];
    ghost[2] = skin[2];
    ghost[3] = skin[3];
    ghost[4] = skin[4];

    // Testing with only momentum at 0.0
    ghost[1] = 0.0;
    ghost[2] = 0.0;
    ghost[3] = 0.0;

    app->in_i_flux = skin[1];
    // printf("Ion in-flux w/ difference: %7.7f\n", skin[1]);
}

//Sheath Ion
void evalIonUpperBC(const struct gkyl_wv_eqn *eqn, double t, int nc, const double *skin, double *GKYL_RESTRICT ghost, void *ctx)
{
    struct sheath_ctx *app = ctx;

    double rhoi_BC =  app->ni * app->mi * 1e-19;
    double Ei_BC = app->ni * app->Ti * 1e-19 / (app->gas_Gamma - 1.0);

    double ghost_mom = ghost[1];
    
    //allowing ions to flow out
    double ui_norm = get_interp_val(&app->data_ui, app->Lx, app->Lx); // find the ui in the rob file 
    double u_local = ui_norm * app->ub; // multiply by the bohm speed to normilize  
    double momi_x = rhoi_BC * u_local; 
    // Set ion upper boundary mass density.
    ghost[0] = rhoi_BC;
    // Set ion upper boundary momentum density.
    ghost[1] = 0.0;
    ghost[2] = 0.0;
    ghost[3] = 0.0;
    // Set ion upper boundary total energy density.
    ghost[4] = Ei_BC;

    app->out_i_flux = skin[1];
    // printf("Ion out-flux w/ difference: %7.7f\n", app->out_i_flux);

    // printf("Value of t (for testing purposes): %7.7f\n", t);
}

void evalElcSourceFunc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx) {
    struct sheath_ctx *app = ctx;
    double x = xn[0];
    int idx = (int)floor((x / app->Lx) * app->Nx);
    double L = app->source_domain_length;
    double nNeut = 1.0;
    double local_ne = 0.0;
    double stangeby = ((2.0 * app->ub) / app->Lx) * (app->pi / 2.0 - 1.0);
    double injection_Te = app->Te / (app->gas_Gamma - 1.0);

    if (app->elc_sheath.data[idx] > 0.0) {
        local_ne = app->elc_sheath.data[idx] / app->me;
    }

    // printf("Value of elc_fout: %7.7f\n", local_ne);
    /* if ((x < L - app->transition_length) && (t > 200.0)) {
        // fout[0] = nNeut * local_ne * stangeby;
        fout[0] = local_ne;
        fout[1] = injection_Te;
    } else if ((x < L) && (t > 200.0)) {
        // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
        double multiplier = (L - x - (app->Lx / app->Nx)) / app->transition_length;
        fout[0] = multiplier * local_ne;
        multiplier = 1.0; // need to remove later
        fout[1] = multiplier * injection_Te;
    } else {
        fout[0] = 0.0;
        fout[1] = injection_Te; // need to change back later
    } */

    double multiplier = 0.0;
    if ((x < L - app->transition_length) && (t > 200.0)) {
        multiplier = 1.0;
    } else if ((x < L) && (t > 200.0)) {
        // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
        multiplier = (L - x) / app->transition_length;
    }

    fout[0] = local_ne;
    fout[1] = injection_Te;
    fout[2] = multiplier;
}

void evalIonSourceFunc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx) {
    struct sheath_ctx *app = ctx;
    double x = xn[0];
    int idx = (int)floor((x / app->Lx) * app->Nx);
    double L = app->source_domain_length;
    double nNeut = 1.0;
    double local_ne = 0.0;
    double stangeby = ((2.0 * app->ub) / app->Lx) * (app->pi / 2.0 - 1.0);
    double injection_Ti = app->Ti / (app->gas_Gamma - 1.0);

    if (app->elc_sheath.data[idx] > 0.0) {
        local_ne = app->elc_sheath.data[idx] / app->me;
    }

    /* if ((x < L - app->transition_length) && (t > 200.0)) {
        // fout[0] = nNeut * local_ne * stangeby;
        fout[0] = local_ne;
        fout[1] = injection_Ti;
    } else if ((x < L) && (t > 200.0)) {
        // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
        double multiplier = (L - x - (app->Lx / app->Nx)) / app->transition_length;
        fout[0] = multiplier * local_ne;
        multiplier = 1.0; // need to remove later
        fout[1] = multiplier * injection_Ti;
    } else {
        fout[0] = 0.0;
        fout[1] = injection_Ti; // need to change back later
    } */

    double multiplier = 0.0;
    if ((x < L - app->transition_length) && (t > 200.0)) {
        multiplier = 1.0;
    } else if ((x < L) && (t > 200.0)) {
        // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
        multiplier = (L - x) / app->transition_length;
    }

    fout[0] = local_ne;
    fout[1] = injection_Ti;
    fout[2] = multiplier;
}

void write_data(struct gkyl_tm_trigger *iot, gkyl_moment_app *app, double t_curr, bool force_write)
{
    if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write)
    {
        int frame = iot->curr - 1;
        if (force_write)
        {
            frame = iot->curr;
        }

        gkyl_moment_app_write(app, t_curr, frame);
        gkyl_moment_app_write_field_energy(app);
        gkyl_moment_app_write_integrated_mom(app);
    }
}

void calc_field_energy(struct gkyl_tm_trigger *fet, gkyl_moment_app *app, double t_curr, bool force_calc)
{
    if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc)
    {
        gkyl_moment_app_calc_field_energy(app, t_curr);
    }
}

void calc_integrated_mom(struct gkyl_tm_trigger *imt, gkyl_moment_app *app, double t_curr, bool force_calc)
{
    if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc)
    {
        gkyl_moment_app_calc_integrated_mom(app, t_curr);
    }
}

int main(int argc, char **argv)
{
    struct gkyl_app_args app_args = parse_app_args(argc, argv);

    #ifdef GKYL_HAVE_MPI
        if (app_args.use_mpi)
        {
            MPI_Init(&argc, &argv);
        }
    #endif

    if (app_args.trace_mem)
    {
        gkyl_cu_dev_mem_debug_set(true);
        gkyl_mem_debug_set(true);
    }

    struct sheath_ctx ctx = create_ctx(); // Context for initialization functions.

    int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

    // Electron/ion equations.
    struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_Gamma, app_args.use_gpu);
    struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_Gamma, app_args.use_gpu);

    struct gkyl_moment_species elc = {
        .name = "elc",
        .charge = ctx.qe,
        .mass = ctx.me,
        .equation = elc_euler,

        .init = evalElcInit,
        .ctx = &ctx,

        .ionization_source_func = evalElcSourceFunc,
        .ionization_source_ctx = &ctx,
        // .ionization_source_set_only_once = true,

        .bcx = {GKYL_SPECIES_FUNC, GKYL_SPECIES_FUNC},
        .bcx_func = {evalElcLowerBC, evalElcUpperBC},
    };
    
    // Ions
    struct gkyl_moment_species ion = {
        .name = "ion",
        .charge = ctx.qi,
        .mass = ctx.mi,
        .equation = ion_euler,
        
        .init = evalIonInit,
        .ctx = &ctx,

        .ionization_source_func = evalIonSourceFunc,
        .ionization_source_ctx = &ctx,
        // .ionization_source_set_only_once = true,

        .bcx = {GKYL_SPECIES_FUNC, GKYL_SPECIES_FUNC},
        .bcx_func = {evalIonLowerBC, evalIonUpperBC},
    };

    // Field.
    struct gkyl_moment_field field = {
        .epsilon0 = ctx.ep0,
        .mu0 = ctx.mu0,

        .is_static = false,
        .init = evalFieldInit,
        .ctx = &ctx,

        .bcx = {GKYL_FIELD_COPY, GKYL_FIELD_PEC_WALL},
    };

    int nrank = 1; // Number of processes in simulation.
    #ifdef GKYL_HAVE_MPI
        if (app_args.use_mpi)
        {
            MPI_Comm_size(MPI_COMM_WORLD, &nrank);
        }
    #endif

    // Create global range.
    int cells[] = {NX};
    int dim = sizeof(cells) / sizeof(cells[0]);
    struct gkyl_range global_r;
    gkyl_create_global_range(dim, cells, &global_r);

    // Create decomposition.
    int cuts[dim];
    #ifdef GKYL_HAVE_MPI
    for (int d = 0; d < dim; d++)
    {
        if (app_args.use_mpi)
        {
            cuts[d] = app_args.cuts[d];
        }
        else
        {
            cuts[d] = 1;
        }
    }
    #else
        for (int d = 0; d < dim; d++)
        {
            cuts[d] = 1;
        }
    #endif

    struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(dim, cuts, &global_r);

    // Construct communicator for use in app.
    struct gkyl_comm *comm;
    #ifdef GKYL_HAVE_MPI
    if (app_args.use_mpi)
    {
        comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){
            .mpi_comm = MPI_COMM_WORLD,
            .decomp = decomp});
    }
    else
    {
        comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
            .decomp = decomp,
            .use_gpu = app_args.use_gpu});
    }
    #else
        comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
            .decomp = decomp,
            .use_gpu = app_args.use_gpu});
    #endif

    int my_rank;
    gkyl_comm_get_rank(comm, &my_rank);
    int comm_size;
    gkyl_comm_get_size(comm, &comm_size);

    int ncuts = 1;
    for (int d = 0; d < dim; d++)
    {
        ncuts *= cuts[d];
    }

    if (ncuts != comm_size)
    {
        if (my_rank == 0)
        {
            fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
        }
        goto mpifinalize;
    }

    // Moment app.
    struct gkyl_moment app_inp = {
        .name = "5m_sheath_60000T_URS_No_Sourse_Instability_dynsrc",

        .ndim = 1,
        .lower = {0.0},
        .upper = {ctx.Lx},
        .cells = {NX},

        .num_periodic_dir = 0,
        .periodic_dirs = {},
        .cfl_frac = ctx.cfl_frac,

        .num_species = 2,
        .species = {elc, ion},

        .has_collision = ctx.has_collision,
        .nu_base = {
            {0.0, ctx.nu_base_ei},
            {ctx.nu_base_ei, 0.0}},

        .field = field,

        .parallelism = {
            .use_gpu = app_args.use_gpu,
            .cuts = {app_args.cuts[0]},
            .comm = comm,
        },
    };

    // Create app object.
    gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

    // Initial and final simulation times.
    double t_curr = 0.0, t_end = ctx.t_end;

    // Initialize simulation.
    int frame_curr = 0;
    if (app_args.is_restart)
    {
        struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

        if (status.io_status != GKYL_ARRAY_RIO_SUCCESS)
        {
            gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
            goto freeresources;
        }

        frame_curr = status.frame;
        t_curr = status.stime;

        gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
        gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
    }
    else
    {
        gkyl_moment_app_apply_ic(app, t_curr);
    }

    // Create trigger for field energy.
    int field_energy_calcs = ctx.field_energy_calcs;
    struct gkyl_tm_trigger fe_trig = {.dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr};

    calc_field_energy(&fe_trig, app, t_curr, false);

    // Create trigger for integrated moments.
    int integrated_mom_calcs = ctx.integrated_mom_calcs;
    struct gkyl_tm_trigger im_trig = {.dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr};

    calc_integrated_mom(&im_trig, app, t_curr, false);

    // Create trigger for IO.
    int num_frames = ctx.num_frames;
    struct gkyl_tm_trigger io_trig = {.dt = t_end / num_frames};

    write_data(&io_trig, app, t_curr, false);

    // Compute estimate of maximum stable time-step.
    double dt = gkyl_moment_app_max_dt(app);

    // Initialize small time-step check.
    double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
    int num_failures = 0, num_failures_max = ctx.num_failures_max;

    long step = 1;

    // Initializing and opening output files
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    FILE *fptr4;
    FILE *fptr5;
    FILE *fptr6;
    FILE *fptr7;
    FILE *fptr8;
    FILE *fptr9;
    fptr1 = fopen("sheath_boundary_ub.txt", "w");
    if (fptr1 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr2 = fopen("elc_particle_count.txt", "w");
    if (fptr2 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr3 = fopen("ion_particle_count.txt", "w");
    if (fptr3 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr4 = fopen("elc_flux_count.txt", "w");
    if (fptr4 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr5 = fopen("ion_flux_count.txt", "w");
    if (fptr5 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr6 = fopen("number_density_benchmark.txt", "w");
    if (fptr6 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr7 = fopen("field_potential_benchmark.txt", "w");
    if (fptr7 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr8 = fopen("ion_velocity_benchmark.txt", "w");
    if (fptr8 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }
    fptr9 = fopen("elc_temp.txt", "w");
    if (fptr9 == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }

    // Defining some key measurement variables
    double elc_tot_flux_count = 0.0;
    double ion_tot_flux_count = 0.0;
    double elc_flux_count = 0.0;
    double ion_flux_count = 0.0;
    double init_elc_count = 0.0;
    double final_elc_count = 0.0;
    double init_ion_count = 0.0;
    double final_ion_count = 0.0;

    double results_spacer = 0.0; // For making sure the benchmarking files only have 100 time-frames

    while ((t_curr < t_end) && (step <= app_args.num_steps))
    {
        gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
        struct gkyl_update_status status = gkyl_moment_update(app, dt);
        gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

        if (!status.success)
        {
            gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
            break;
        }

        t_curr += status.dt_actual;
        dt = status.dt_suggested;

        calc_field_energy(&fe_trig, app, t_curr, false);
        calc_integrated_mom(&im_trig, app, t_curr, false);
        write_data(&io_trig, app, t_curr, false);

        if (dt_init < 0.0)
        {
            dt_init = status.dt_actual;
        }
        else if (status.dt_actual < dt_failure_tol * dt_init)
        {
            num_failures += 1;

            gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
            gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
            gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
            if (num_failures >= num_failures_max)
            {
                gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
                gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

                calc_field_energy(&fe_trig, app, t_curr, true);
                calc_integrated_mom(&im_trig, app, t_curr, true);
                write_data(&io_trig, app, t_curr, true);
                break;
            }
        }
        else
        {
            num_failures = 0;
        }

        // Redefining fluxes
        const double *f_elc_a = gkyl_array_fetch(app->species[0].f[0], 0);
        const double *f_ion_a = gkyl_array_fetch(app->species[1].f[0], 0);
        const double *f_elc_b = gkyl_array_fetch(app->species[0].f[0], NX - 1);
        const double *f_ion_b = gkyl_array_fetch(app->species[1].f[0], NX - 1);
        ctx.out_e_flux = f_elc_b[1] - f_elc_a[1];
        ctx.out_i_flux = f_ion_b[1] - f_ion_a[1];

        // Update scheme that bypasses nT source
        /*
        for (int idx = 0; idx < NX; idx++) {
            double *fw_elc = gkyl_array_fetch(app->species[0].f[0], idx);
            double *fw_ion = gkyl_array_fetch(app->species[1].f[0], idx);

            double x = idx * (ctx.Lx / NX);
            double t = t_curr;
            double deltaT = status.dt_actual;
            double L = ctx.source_domain_length;
            double local_ne = 0.0;
            double local_me_old = fw_elc[0];
            double local_mi_old = fw_ion[0];
            double injection_Te = ctx.Te;
            double injection_Ti = ctx.Ti;
            double prev_Te = (fw_elc[4] - 0.5 * fw_elc[1] * fw_elc[1] / fw_elc[0]) / (fw_elc[0] / ctx.me);
            double prev_Ti = (fw_ion[4] - 0.5 * fw_ion[1] * fw_ion[1] / fw_ion[0]) / (fw_ion[0] / ctx.mi);

            if (ctx.elc_sheath.data[idx] > 0.0) {
                local_ne = ctx.elc_sheath.data[idx] / ctx.me;
            }
            if ((x < L - ctx.transition_length) && (t > 200.0)) {
                fw_elc[0] += deltaT * local_ne * ctx.me;
                fw_ion[0] += deltaT * local_ne * ctx.mi;
            } else if ((x < L) && (t > 200.0)) {
                // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
                double multiplier = (L - x - (ctx.Lx / NX)) / ctx.transition_length;
                // fw_elc[0] += deltaT * multiplier * ctx.elc_sheath.data[(int)floor(((L - 4.0 - ctx.Lx / NX) / ctx.Lx) * NX)];
                // fw_ion[0] += deltaT * multiplier * ctx.elc_sheath.data[(int)floor(((L - 4.0 - ctx.Lx / NX) / ctx.Lx) * NX)] * (ctx.mi / ctx.me);
                fw_elc[0] += deltaT * multiplier * local_ne * ctx.me;
                fw_ion[0] += deltaT * multiplier * local_ne * ctx.mi;
            }

            if (t > 200.0) {
                // No changes to f_elc[1] or f_ion[0] since ionized neutrals are added at rest
                double added_me = fw_elc[0] - local_me_old;
                double added_mi = fw_ion[0] - local_mi_old;

                // In order to reincorporate adiabatic constant based on its presence in init functions
                injection_Te = injection_Te / (ctx.gas_Gamma - 1.0);
                injection_Ti = injection_Ti / (ctx.gas_Gamma - 1.0);

                // While the kinetic energy term would decrease due to the added particles and constant momentum, the energy would be conserved as increased temperature.
                // fw_elc[4] = ((local_me_old * prev_Te + added_me * injection_Te) / ctx.me) + (0.5 * fw_elc[1] * fw_elc[1] / local_me_old);
                // fw_ion[4] = ((local_mi_old * prev_Ti + added_mi * injection_Ti) / ctx.mi) + (0.5 * fw_ion[1] * fw_ion[1] / local_mi_old);

                // Alternative: sets temperature to constant for all fluids. May not be so good for smaller domain size, or may require being applied only across source term domain.
                if (x < L - ctx.transition_length) {
                    fw_elc[4] = ((local_me_old * injection_Te + added_me * injection_Te) / ctx.me) + (0.5 * fw_elc[1] * fw_elc[1] / fw_elc[0]);
                    fw_ion[4] = ((local_mi_old * injection_Ti + added_mi * injection_Ti) / ctx.mi) + (0.5 * fw_ion[1] * fw_ion[1] / fw_ion[0]);
                } else if (x < L) {
                    // Transition region in order to prevent instability formation at source term domain edge. Tradeoff is particle flux accuracy.
                    double multiplier = (L - x - (ctx.Lx / NX)) / ctx.transition_length;
                    fw_elc[4] += multiplier * (((local_me_old * injection_Te + added_me * injection_Te) / ctx.me) + (0.5 * fw_elc[1] * fw_elc[1] / fw_elc[0]) - fw_elc[4]);
                    fw_ion[4] += multiplier * (((local_mi_old * injection_Ti + added_mi * injection_Ti) / ctx.mi) + (0.5 * fw_ion[1] * fw_ion[1] / fw_ion[0]) - fw_ion[4]);
                }
            }
        }
        */

        // For determining boundary between pre-sheath and sheath, as well as writing total particle count to file
        double ion_fvel_diff = -1.0; // dummy value, replaced almost immediately
        double ion_min_fvel_x = 0.0;
        double elc_particle_count = 0.0;
        double ion_particle_count = 0.0;
        double halfd_elc_particle_count = 0.0;
        double halfd_ion_particle_count = 0.0;
        double tot_potential = 0.0;
        int source_domain_length_NX = ((int)floor(ctx.source_domain_length * NX / ctx.Lx));

        // Since we replaced elc flux with ion flux earlier
        // elc_flux_count = (ctx.out_e_flux - ctx.in_e_flux) * dt;
        elc_flux_count = (ctx.out_i_flux - ctx.in_i_flux) * dt * (ctx.me / ctx.mi);

        // For testing to see if instabilities carry over
        ion_flux_count = (ctx.out_i_flux - ctx.in_i_flux) * dt;
        // ion_flux_count = (ctx.out_e_flux - ctx.in_e_flux) * dt * (ctx.mi / ctx.me);
        elc_tot_flux_count += elc_flux_count;
        ion_tot_flux_count += ion_flux_count;
        for (int sel = 0; sel < NX; sel++) {
            // NX is number of cells in x-direction, i.e. number of cells in system since this is 1D
            int sidx[] = { 0, app->ndim };
            const double *f_elc = gkyl_array_fetch(app->species[0].f[sidx[0]], sel);
            const double *f_ion = gkyl_array_fetch(app->species[1].f[sidx[0]], sel);
            double elc_fvel = f_elc[1] / f_elc[0];
            double ion_fvel = f_ion[1] / f_ion[0];
            // printf("Value of elc_fvel (for testing purposes): %7.7f\n", elc_fvel);
            // printf("Value of ion_fvel (for testing purposes): %7.7f\n", ion_fvel);
            if (sel <= 0) {
                ion_fvel_diff = fabs(ion_fvel - ctx.ub);
            }
            if (fabs(ion_fvel - ctx.ub) < ion_fvel_diff) {
                ion_fvel_diff = fabs(ion_fvel - ctx.ub);
                ion_min_fvel_x = sel * (ctx.Lx / NX);
            }

            elc_particle_count += f_elc[0] * (ctx.Lx / NX); // This and ion particle count are basically just the integral of their respective number densities over the domain, expressed here as riemann sums
            // Okay so mass densities, not number densities. Which means the nomenclature doesn't really make sense so later you should probably divide both of these by their respective particle masses and then modify the python sheath grapher code
            ion_particle_count += f_ion[0] * (ctx.Lx / NX);

            // Attempting to access field data
            const double *f_field = gkyl_array_fetch(app->field.f[sidx[0]], sel);
            tot_potential += f_field[0];
            if (sel < source_domain_length_NX - (int)(ctx.transition_length * NX / ctx.Lx)) {
                halfd_elc_particle_count += f_elc[0] * (ctx.Lx / NX);
                halfd_ion_particle_count += f_ion[0] * (ctx.Lx / NX);
            } else if (sel < source_domain_length_NX) {
                double multiplier = (source_domain_length_NX - sel - 1) / (ctx.transition_length * NX / ctx.Lx);
                /* int transition_sel = source_domain_length_NX - (int)(4.0 * NX / ctx.Lx) - 1;
                const double *f_elc_transition = gkyl_array_fetch(app->species[0].f[sidx[0]], transition_sel);
                const double *f_ion_transition = gkyl_array_fetch(app->species[1].f[sidx[0]], transition_sel); */
                halfd_elc_particle_count += f_elc[0] * multiplier * (ctx.Lx / NX);
                halfd_ion_particle_count += f_ion[0] * multiplier * (ctx.Lx / NX);
            }
        }

        // Creating runtime sheath profiles for particle injection
        double elc_sheath_array[NX];
        double ion_sheath_array[NX];
        for (int sel = 0; sel < source_domain_length_NX; sel++) {
            int sidx[] = { 0, app->ndim };
            const double *f_elc = gkyl_array_fetch(app->species[0].f[sidx[0]], sel);
            const double *f_ion = gkyl_array_fetch(app->species[1].f[sidx[0]], sel);
            double collisionality_operator = ((ctx.out_i_flux - ctx.in_i_flux) * (ctx.me / ctx.mi) / halfd_elc_particle_count);
            // double collisionality_operator = 0.00018; // Setting to constant to see what happens
            if (collisionality_operator > 0.0) {
                elc_sheath_array[sel] = f_elc[0] * collisionality_operator;
            } else {
                elc_sheath_array[sel] = 0.0;
            }
            ion_sheath_array[sel] = f_ion[0];
            ctx.elc_sheath.data[sel] = elc_sheath_array[sel];
            ctx.ion_sheath.data[sel] = ion_sheath_array[sel];
        }
        for (int sel = source_domain_length_NX; sel < NX; sel++) {
            ctx.elc_sheath.data[sel] = 0.0;
            ctx.ion_sheath.data[sel] = 0.0;
        }
        // const double *f_elc = gkyl_array_fetch(app->species[0].f[0], 1);
        // const double *f_ion = gkyl_array_fetch(app->species[1].f[0], 1);
        // printf("Value of ion_sheath.data[1]: %7.7f\n", ctx.ion_sheath.data[1] / f_ion[0]);
        // printf("Value of elc_sheath.data[1]: %7.7f\n", ctx.elc_sheath.data[1] / f_elc[0]);

        // Reporting values to console and saving some to ctx
        ctx.presheath_length = ion_min_fvel_x;
        printf("Value of ion_min_fvel_x: %7.7f\n", ion_min_fvel_x);
        printf("Value of elc_particle_count: %7.7f\n", elc_particle_count);
        printf("Value of ion_particle_count: %7.7f\n", ion_particle_count);
        printf("Value of halfd_elc_particle_count: %7.7f\n", halfd_elc_particle_count);
        printf("Value of halfd_ion_particle_count: %7.7f\n", halfd_ion_particle_count);
        printf("Value of elc_particle_count difference: %7.7f\n", (elc_particle_count - ctx.e_count) * (1.0 / ctx.me));
        printf("Value of ion_particle_count difference: %7.7f\n", (ion_particle_count - ctx.i_count) * (1.0 / ctx.mi));
        printf("Value of elc number flux at right wall: %7.7f\n", (ctx.out_e_flux - ctx.in_e_flux) * dt * (1.0 / ctx.me));
        printf("Value of ion number flux at right wall: %7.7f\n", (ctx.out_i_flux - ctx.in_i_flux) * dt * (1.0 / ctx.mi));
        int sidx[] = { 0, app->ndim };
        const double *f_elc = gkyl_array_fetch(app->species[0].f[sidx[0]], ion_min_fvel_x * (int)(NX / ctx.Lx));
        const double *f_ion = gkyl_array_fetch(app->species[1].f[sidx[0]], ion_min_fvel_x * (int)(NX / ctx.Lx));
        printf("Value of elc number flux at sheath boundary: %7.7f\n", f_elc[0] * ctx.ub * dt * (1.0 / ctx.me));
        printf("Value of ion number flux at sheath boundary: %7.7f\n", f_ion[0] * ctx.ub * dt * (1.0 / ctx.mi));
        printf("Value of ionization source collisionality operator: %7.7f\n", ((ctx.out_i_flux - ctx.in_i_flux) * (ctx.me / ctx.mi) / halfd_elc_particle_count));
        ctx.e_count = elc_particle_count;
        ctx.i_count = ion_particle_count;

        // writing values to files
        fprintf(fptr1, "%7.7f\n", ion_min_fvel_x);
        fprintf(fptr2, "%7.7f\n", elc_particle_count);
        fprintf(fptr3, "%7.7f\n", ion_particle_count);
        fprintf(fptr4, "%7.7f\n", elc_particle_count + elc_tot_flux_count);
        fprintf(fptr5, "%7.7f\n", ion_particle_count + ion_tot_flux_count);
        // This bit here is for getting the initial and final particle counts in system for printing to console
        if (step <= 1) {
            init_elc_count = elc_particle_count;
            init_ion_count = ion_particle_count;
        }
        if ((t_curr >= t_end) || (step >= app_args.num_steps)) {
            final_elc_count = elc_particle_count;
            final_ion_count = ion_particle_count;
        }

        // For benchmarking files
        if (t_curr >= results_spacer) {
            int NX_multiplier = (int)floor(NX / ctx.Lx);
            for (int i = (NX - 28 * NX_multiplier); i < NX; i++) {
                int sidx[] = { 0, app->ndim };
                const double *f_elc = gkyl_array_fetch(app->species[0].f[sidx[0]], i);
                const double *f_ion = gkyl_array_fetch(app->species[1].f[sidx[0]], i);
                double local_potential = tot_potential;
                for (int j = 0; j <= i; j++) {
                    const double *f_field = gkyl_array_fetch(app->field.f[sidx[0]], j);
                    local_potential -= f_field[0];
                }
                const double *f_field = gkyl_array_fetch(app->field.f[sidx[0]], i);
                double local_field = f_field[0];
                double ion_fvel_relative = (f_ion[1] / f_ion[0]) / ctx.ub;
                double num_density_e = f_elc[0] / ctx.me;
                double num_density_i = f_ion[0] / ctx.mi;
                double num_density_diff = num_density_i - num_density_e;

                fprintf(fptr6, "%7.7f, ", num_density_e);
                fprintf(fptr6, "%7.7f, ", num_density_i);
                fprintf(fptr6, "%7.7f\n", num_density_diff);
                fprintf(fptr7, "%7.7f, ", local_potential);
                fprintf(fptr7, "%7.7f\n", local_field);
                fprintf(fptr8, "%7.7f\n", ion_fvel_relative);
            }
            for (int i = 0; i < NX; i++) {
                int sidx[] = { 0, app->ndim };
                const double *f_elc = gkyl_array_fetch(app->species[0].f[sidx[0]], i);
                double elc_temp = (f_elc[4] - 0.5 * f_elc[1] * f_elc[1] / f_elc[0]) / (f_elc[0] / ctx.me);
                fprintf(fptr9, "%7.7f\n", elc_temp);
            }
            fprintf(fptr6, "\n");
            fprintf(fptr7, "\n");
            fprintf(fptr8, "\n");
            fprintf(fptr9, "\n");
            results_spacer += 200.0;
            printf("results_spacer updated to %7.7f\n", results_spacer);
        }

        step += 1;

        // Clear terminal
        /* printf("\e[1;1H\e[2J");
        fflush(stdout); */
    }

    // Closing files
    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    fclose(fptr4);
    fclose(fptr5);
    fclose(fptr6);
    fclose(fptr7);
    fclose(fptr8);
    fclose(fptr9);

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);
    gkyl_moment_app_stat_write(app);

    struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

    gkyl_moment_app_cout(app, stdout, "\n");
    gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
    gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
    gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
    gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
    gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
    gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

    // Printing out lost particles from flux for future comparison
    printf("Total lost electron flux %7.7f\n", elc_tot_flux_count);
    printf("Total lost ion flux %7.7f\n", ion_tot_flux_count);
    double lost_elc_count = -(final_elc_count - init_elc_count);
    double lost_ion_count = -(final_ion_count - init_ion_count);
    printf("Total lost electron mass %7.7f\n", lost_elc_count);
    printf("Total lost ion mass %7.7f\n", lost_ion_count);

    // Free resources after simulation completion.
    freeresources:
    gkyl_wv_eqn_release(elc_euler);
    gkyl_wv_eqn_release(ion_euler);
    gkyl_rect_decomp_release(decomp);
    gkyl_comm_release(comm);
    gkyl_moment_app_release(app);

    mpifinalize:
    #ifdef GKYL_HAVE_MPI
        if (app_args.use_mpi)
        {
            MPI_Finalize();
        }
    #endif

    return 0;
}
