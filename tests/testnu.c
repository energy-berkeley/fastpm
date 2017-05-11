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

#include "../libfastpm/pmpfft.h"

//#include <fastpm/cosmology.h>
//We would like to couple neutrinos into the dark matter evolution.
//Stage 1: read in and record all the overdensity field. ---- done
//Let's do Simpson rule instead. So no interpolation needed. -- to stage 3
 
//Stage 2: interpolate the read overdensity. ---- skipped first.
//Stage 3: integrate out the overdendity --- in progress.
//Stage 4: feed the neutrino overdensity into force calculation

//check units and numbers, k is absolutely not right.
// Side check: time fourier transform of the delta k, frequncy amplitude to look at smothness.
// Change to check with linear theory.

#define hhhh  70
#define m  0.6; //in eV
#define c  3e8; // in m/s
#define chunit  1.68e-4/0.6/0.678*3e8/1000; //used m, hubble and c

typedef struct {
    FastPMSolver * solver;
    FastPMFloat ** tape;
    FastPMFloat ** Nu;
    PM * pm;
    int step;
    int maxsteps;
    double * time_step;
    double *SupCo;
} FastPMRecorder;

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps, double * time_step);

void fastpm_recorder_destroy(FastPMRecorder * recorder);

static void
record_cdm(FastPMSolver * solver, FastPMForceEvent * event, FastPMRecorder * recorder);

static void Del_interp(FastPMRecorder * recorder);

double Sint(double a, void * params){// FastPMSolver * solver){
    FastPMCosmology * cosmo = (FastPMCosmology *) params;
    return 1.0/a/a/a/HubbleEa(a, cosmo)/hhhh;  //How to implement this correctly? HubbleEa(a, fastpm->cosmology)*71.9; //70.0; //Planck15.H((1.-a)/a)
}

double SupCon(double ai,double af, FastPMSolver * solver){
    gsl_integration_workspace * w= gsl_integration_workspace_alloc (100);
    double result, error;
    gsl_function F;

    F.function = &Sint;
    F.params = &solver->cosmology;
    gsl_integration_qags(&F, ai, af, 0, 1e-7, 100, w, &result, &error);
    return result;
}

double kdifs(double k, FastPMRecorder * recorder, int j){
    double difs = recorder->SupCo[j];//SupCon(a,af,solver);
    return k*difs;
}

double CurlyInum(double x){
    return (1.0+0.0168*x*x+0.0407*pow(x,4));
}

double CurlyIden(double x){
    return 1.0+2.1734*x*x+1.6787*pow(x,4.1811)+0.1467*pow(x,8);
}

double CurlyI(double x){
    return CurlyInum(x)/CurlyIden(x);
}

double Inte(double k,FastPMRecorder * recorder,int j){
    double a = recorder->time_step[j];
    double x = kdifs(k,recorder,j)*chunit;
    return CurlyI(x)*kdifs(k,recorder,j)/k/HubbleEa(a, recorder->solver->cosmology)/hhhh/a/a;// Planck15.H((1.-a)/a).value/a**2
}

double simpson(double * data, double da, int i);

double DelNu(double * Del,double k, FastPMRecorder * recorder){
    double data[recorder->maxsteps];//FIXME hard coded array size,fixed!
    int j = 0;
    for(;j<=recorder->step;++j){
        data[j] = Inte(k,recorder,j)*Del[j];/*af = a[9]? Checked python, ok.*/
    }
    double result = simpson(data, recorder->time_step[1]-recorder->time_step[0], recorder->step);
    return result;
}



int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    double time_step[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, .9, 1.0};
//    double time_step[] = {0.0  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
//        0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
//        0.9 ,  0.95,  1.0};

    FastPMConfig * config = & (FastPMConfig) {
        .nc = 8,
        .boxsize = 8.,
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
        .Norm = 3000.0, /* FIXME: this is not any particular sigma8. */
        .hubble_param = hhhh/100.0,
        .omegam = 0.260,
        .omegab = 0.044,
    };
    fastpm_ic_fill_gaussiank(solver->basepm, rho_init_ktruth, 2004, FASTPM_DELTAK_GADGET);
    fastpm_ic_induce_correlation(solver->basepm, rho_init_ktruth, (fastpm_fkfunc)fastpm_utils_powerspec_eh, &eh);

    write_complex(solver->basepm, rho_init_ktruth, "Truth", "Delta_k", 1);
    fastpm_solver_setup_ic(solver, rho_init_ktruth);

    fastpm_recorder_init(recorder, solver, 10, time_step);

    fastpm_solver_evolve(solver, time_step, sizeof(time_step) / sizeof(time_step[0]));


    printf("I am running\n");
    

    pm_free(solver->basepm, rho_init_ktruth);
    fastpm_recorder_destroy(recorder);

    fastpm_solver_destroy(solver);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}

void fastpm_recorder_init(FastPMRecorder * recorder, FastPMSolver * solver, int maxsteps, double * time_step)
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
    recorder->time_step = time_step;
    double SuCo[maxsteps];
    int l; 
    for(l =0; l <= (maxsteps-1); ++l){
        SuCo[l] = SupCon(time_step[l],time_step[maxsteps-1], solver);
        }
    recorder->SupCo = SuCo;
    int mo=0;
    for(;mo<=maxsteps-1;++mo)printf("SupCon is %g\n", SuCo[mo]);
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
    sprintf(buf, "Nu%0.04f.dat", event->a_f);
    printf("The step %g\n", event->a_f);
    FastPMFloat * dst = recorder->tape[recorder->step];
    PM * pm = solver->pm;
    pm_assign(pm, event->delta_k, dst);
    PMKIter kiter;
    for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {

        int d = 0;
        ptrdiff_t ind = kiter.ind;

            ptrdiff_t kk = 0.;
            for(d = 0; d < 3; d++) {
                double ik = kiter.iabs[d];
                if(ik > pm->Nmesh[d] / 2) ik -= pm->Nmesh[d];
                kk += ik * ik;
            }

            double scark = sqrt(kk) * 2 * M_PI / pm_boxsize(pm)[0];
            // Start integrating
            int j = 0;
            double realDel[recorder->maxsteps], ImaDel[recorder->maxsteps];
            double realNu[recorder->maxsteps], ImaNu[recorder->maxsteps];
            for(j=0;j<=recorder->step;++j){
                realDel[j] = recorder->tape[j][ind + 0];
                ImaDel[j] = recorder->tape[j][ind + 1];
                recorder->Nu[j][ind + 0] = DelNu(realDel,scark,recorder);
                recorder->Nu[j][ind + 1] = DelNu(ImaDel,scark,recorder);
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


double simpson(double * data, double da, int i){

    double inte = 0;
    int j = 0;
    for(j =0;j<=i;++j){
        if(!(i%2)){
            if(j == 0) inte = da/3*data[0];
            else if (j>0 && j%2 && j<i) inte += 4*da/3*data[j];
            else if (j>0 && !(j%2) && j<i) inte += 2*da/3*data[j];
            else if (j == i) inte += da/3*data[j];
        }
        else{
            if(j == 0) inte = da/3*data[0];
            else if (j>0 && j%2 && j<i) inte += 2*da/3*data[j];
            else if (j>0 && !(j%2) && j<i) inte += 4*da/3*data[j];
            else if (j == i) inte += da/3*data[j];
        }
    }
    return inte;
}

