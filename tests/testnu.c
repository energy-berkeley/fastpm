#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

//We would like to couple neutrinos into the dark matter evolution.
//Stage 1: read in and record all the overdensity field. ---- done
//Let's do Simpson rule instead. So no interpolation needed. -- to stage 3
 
//Stage 2: interpolate the read overdensity. ---- skipped first.
//Stage 3: integrate out the overdendity --- in progress.
//Stage 4: feed the neutrino overdensity into force calculation


// Side check: time fourier transform of the delta k, frequncy amplitude to look at smothness.
// Change to check with linear theory.

typedef struct {
    FastPMSolver * solver;
    FastPMFloat ** tape;
    FastPMFloat ** Nu;
    PM * pm;
    int step;
    int maxsteps;
} FastPMRecorder;

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps);

void fastpm_recorder_destroy(FastPMRecorder * recorder);

static void
record_cdm(FastPMSolver * solver, FastPMForceEvent * event, FastPMRecorder * recorder);

static void Del_interp(FastPMRecorder * recorder);

double Sint(double a){
    return 1.0/pow(a,3)/20.0;  //How to implement this correctly? HubbleEa(a, fastpm->cosmology); //70.0; //Planck15.H((1.-a)/a)
}

double SupCon(double ai,double af){
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (1000);
    double result, error;
    gsl_function F;

    F.function = &Sint;
    gsl_integration_qags(&F, ai, af, 1e-7, 1000, w, &result, &error);
    return result;
}

double kdifs(double *k,double ai,double af,double a){
    return kdif
kdif[j] = k[j]*SupCon(a[j],af);
}

double CurlyInum(double *x){
    return (1.+0.0168*x^2+0.0407*x^4);
}

double CurlyIden(double *x){
    return 1.+2.1734*x^2+1.6787*x^4.1811+0.1467*x^8;
}

double CurlyI(double *x){
    int j = 0;
    double I[] = x;
    for(;x[j];++j) I[j] = CurlyInum(x[j])/CurlyIden(x[j]);
    return I;
}


double simpson(double ** data, double da, int i, int ind);


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 256,
        .boxsize = 256.,
        .alloc_factor = 2.0,
        .omega_m = 0.292,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = 2},
            {.a_start = -1, .pm_nc_factor = 0},
        },
        .FORCE_TYPE = FASTPM_FORCE_FASTPM,
        .nLPT = 2.5,
        .K_LINEAR = 0.04,
        .SAVE_Q = 1,
    };


    FastPMSolver solver[1];
    FastPMRecorder recorder[1];
    fastpm_solver_init(solver, config, comm);

    FastPMFloat * rho_init_ktruth = pm_alloc(solver->basepm);

    fastpm_solver_add_event_handler(solver, FASTPM_EVENT_FORCE,
            FASTPM_EVENT_STAGE_BEFORE,
            (FastPMEventHandlerFunction) record_cdm,
            recorder);


    /* First establish the truth by 2lpt -- this will be replaced with PM. */
    struct fastpm_powerspec_eh_params eh = {
        .Norm = 10000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = 0.7,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    write_complex(solver->basepm, rho_init_ktruth, "Truth", "Delta_k", 1);
    double time_step[] = {0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
//    double time_step[] = {0.0  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
//        0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
//        0.9 ,  0.95,  1.0};
    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_recorder_init(recorder, solver, 10);

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));


    printf("I am running\n");
    

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_recorder_destroy(recorder);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps)
{
    recorder->maxsteps = maxsteps;
    recorder->solver = solver;
    recorder->tape = calloc(maxsteps, sizeof(FastPMStore));
    recorder->Nu = calloc(maxsteps, sizeof(FastPMStore));
	recorder->step = 0;
    recorder->pm = fastpm_find_pm(solver, 1.0);
    int j =0;
    for(j =0; j <= (maxsteps-1); ++j){
        recorder->tape[j] = pm_alloc(recorder->pm);
        recorder->Nu[j] = pm_alloc(recorder->pm);
        }
}


void fastpm_recorder_destroy(FastPMRecorder * recorder)
{
    int j =0;
    for(j =0; j <= (recorder->maxsteps-1); ++j){
        pm_free(recorder->pm, recorder->tape[j]);
        pm_free(recorder->pm, recorder->Nu[j]);
        }
    free(recorder->tape);
    free(recorder->Nu);
}




static void record_cdm(FastPMSolver * solver, FastPMForceEvent * event, FastPMRecorder * recorder)
{
    char buf[1024];
    sprintf(buf, "Smooth%0.04f.dat", event->a_f);
    printf("The step %g\n", event->a_f);
    FastPMFloat * dst = recorder->tape[recorder->step];
    pm_assign(solver->pm, event->delta_k, dst);
    PMKIter kiter;
    for(pm_kiter_init(solver->pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {

        ptrdiff_t ind = kiter.ind;
        if(kiter.iabs[0] == 1 &&
                kiter.iabs[1] == 1 &&
                kiter.iabs[2] == 1) {
            double real1 = event->delta_k[ind + 0];
            double imag1 = event->delta_k[ind + 1];
            double value1 = real1 * real1 + imag1 * imag1;
            if(recorder->tape[recorder->step]) printf("The tape was initialized!\n");
            double real2 = recorder->tape[recorder->step][ind + 0];
            double imag2 = recorder->tape[recorder->step][ind + 1];
            double value2 = real2 * real2 + imag2 * imag2;
            printf("Delta_k from event = %g + %gi, abs = %g\n", real1, imag1, value1);
            printf("Delta_k from tape = %g + %gi, abs = %g\n", real2, imag2, value2);
            // Start integrating
            int j = 0;
            
            for(j=0;j<=recorder->step;++j){
                recorder->Nu[recorder->step][ind + 0] = simpson(recorder->tape, 0.1, recorder->step, ind);
                recorder->Nu[recorder->step][ind + 1] = simpson(recorder->tape, 0.1, recorder->step, ind+1);
            }
            double real3 = recorder->Nu[recorder->step][ind + 0];
            double imag3 = recorder->Nu[recorder->step][ind + 1];
            double value3 = real3 * real3 + imag3 * imag3;
            printf("Nu = %g + %gi, abs = %g\n", real3, imag3, value3);


        }

    } 
    recorder->step = recorder->step +1;
    write_complex(solver->pm, event->delta_k, buf, "Delta_k", 1);
}

static void Del_interp(FastPMRecorder * recorder)
{
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    const size_t N = 100;
// we need to interpolate all kx ky kz????? 
// N^3 interpolation?
}


double simpson(double ** data, double da, int i, int ind){

    double inte = 0;
    int j = 0;
    for(j =0;j<=i;++j){
        if(!(i%2)){
            if(j == 0) inte = da/3*data[0][ind];
            else if (j>0 && j%2 && j<i) inte += 4*da/3*data[j][ind];
            else if (j>0 && !(j%2) && j<i) inte += 2*da/3*data[j][ind];
            else if (j == i) inte += da/3*data[j][ind];
        }
        else{
            if(j == 0) inte = da/3*data[0][ind];
            else if (j>0 && j%2 && j<i) inte += 2*da/3*data[j][ind];
            else if (j>0 && !(j%2) && j<i) inte += 4*da/3*data[j][ind];
            else if (j == i) inte += da/3*data[j][ind];
        }
    }
    return inte;
}

