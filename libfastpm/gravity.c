#include <string.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/transfer.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>

#include "pmpfft.h"
#include "pmghosts.h"

static void
apply_pot_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int order)
{
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        float ** kklist [3] = {kiter.kk, kiter.kk_finite, kiter.kk_finite2};
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kklist[order][d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - 1 / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] = - from[ind + 0] * (1 / kk_finite);
                to[ind + 1] = - from[ind + 1] * (1 / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

static void
apply_grad_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir, int order)
{
#pragma omp parallel 
    {
        PMKIter kiter;
        ptrdiff_t * Nmesh = pm_nmesh(pm);
        pm_kiter_init(pm, &kiter);
        float ** klist[2] = {kiter.k, kiter.k_finite};
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k_finite;
            k_finite = klist[order][dir][kiter.iabs[dir]];
            ptrdiff_t ind = kiter.ind;
            /* i k[d] */
            /* Watch out the data dependency */
            if(
                kiter.iabs[0] == (Nmesh[0] - kiter.iabs[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - kiter.iabs[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - kiter.iabs[2]) % Nmesh[2]
            ) {
                /* We are at the nyquist and the diff operator shall be zero;
                 * otherwise the force is not real! */
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            } else {
                FastPMFloat tmp = from[ind + 0] * (k_finite);
                to[ind + 0] = - from[ind + 1] * (k_finite);
                to[ind + 1] = tmp;
            }
        }
    }
}

static void
apply_gaussian_dealiasing(PM * pm, FastPMFloat * from, FastPMFloat * to, double N)
{
    /* N is rms in mesh size */
    double r0 = N * pm->BoxSize[0] / pm->Nmesh[0];

#pragma omp parallel 
    {
        PMKIter kiter;
        int d;
        int i;
        double *kernel[3];
        pm_kiter_init(pm, &kiter);
        for(d = 0; d < 3; d ++) {
            kernel[d] = malloc(sizeof(double) * pm->Nmesh[d]);
            for(i = 0; i < pm->Nmesh[d]; i ++) {
                kernel[d][i] = exp(- 0.5 * pow(kiter.k[d][i] * r0, 2));
            }
        }

        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double fac = 1;
            ptrdiff_t ind = kiter.ind;
            for(d = 0; d < 3; d++) {
                fac *= kernel[d][kiter.iabs[d]];
            }
            to[ind + 0] *= fac;
            to[ind + 1] *= fac;
        }
    }
}
static double
gaussian36(double k, double * knq)
{
    double x = k / *knq;
    return exp(- 36 * pow(x, 36));
}

static void
apply_kernel_transfer(FastPMGravity * gravity, PM * pm, FastPMFloat * delta_k, FastPMFloat * canvas, int d)
{
    switch(gravity->KernelType) {
        case FASTPM_KERNEL_EASTWOOD:
            apply_pot_transfer(pm, delta_k, canvas, 0);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 0);
            /* now sharpen for mass assignment */
            /* L1 */
            fastpm_apply_decic_transfer(pm, canvas, canvas);
            /* L2 */
            fastpm_apply_decic_transfer(pm, canvas, canvas);
        break;
        case FASTPM_KERNEL_NAIVE:
            apply_pot_transfer(pm, delta_k, canvas, 0);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 0);
        break;
        case FASTPM_KERNEL_GADGET:
            apply_pot_transfer(pm, delta_k, canvas, 0);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 1);
        break;
        case FASTPM_KERNEL_3_4:
            apply_pot_transfer(pm, delta_k, canvas, 1);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 1);
        break;
        case FASTPM_KERNEL_5_4:
            apply_pot_transfer(pm, delta_k, canvas, 2);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 1);
        break;
        case FASTPM_KERNEL_3_2:
            apply_pot_transfer(pm, delta_k, canvas, 1);
            if(d != 3)
                apply_grad_transfer(pm, canvas, canvas, d, 0);
        break;
        default:
            fastpm_raise(-1, "Wrong kernel type\n");
    }
}
static void
apply_dealiasing_transfer(FastPMGravity * gravity, PM * pm, FastPMFloat * from, FastPMFloat * to)
{
    switch(gravity->DealiasingType) {
        case FASTPM_DEALIASING_TWO_THIRD:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_lowpass_transfer(pm, from, to, 2.0 / 3 * k_nq);
            }
        break;
        case FASTPM_DEALIASING_GAUSSIAN:
            apply_gaussian_dealiasing(pm, from, to, 1.0);
        break;
        case FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN:
            apply_gaussian_dealiasing(pm, from, to, 4.0);
        break;
        case FASTPM_DEALIASING_GAUSSIAN36:
            {
            double k_nq = M_PI / pm->BoxSize[0] * pm->Nmesh[0];
            fastpm_apply_any_transfer(pm, from, to, (fastpm_fkfunc) gaussian36, &k_nq);
            }
        break;
        case FASTPM_DEALIASING_NONE:
        break;
        default:
            fastpm_raise(-1, "wrong dealiasing kernel type");
    }
}

void
fastpm_solver_gravity_calculate(FastPMSolver * fastpm,
    FastPMGravity * gravity,
    PM * pm,
    FastPMStore * p,
    FastPMTransition * trans)
{
    FastPMFloat * delta_k = pm_alloc(pm);
    FastPMPainter reader[1];
    FastPMPainter painter[1];

    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);
    fastpm_painter_init(painter, pm, gravity->PainterType, gravity->PainterSupport);

    /* watch out: boost the density since mesh is finer than grid */
    long long np = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, pm_comm(pm));

    double density_factor = pm->Norm / np;

    CLOCK(ghosts);
    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);
    LEAVE(ghosts);

    FastPMFloat * canvas = pm_alloc(pm);

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * We thus have to boost the density by density_factor.
     *
     * This gives us over density + 1
     *
     * because rhobar = N_g ^3 / V
     * we paint rho V / (B N_g^3) * B = rho / rhobar. The last B is the extra density factor.
     * */
    CLOCK(paint);
    fastpm_paint_local(painter, canvas,
                p, p->np + pgd->nghosts, NULL, 0);
    fastpm_apply_multiply_transfer(pm, canvas, canvas, density_factor);
    LEAVE(paint);

    CLOCK(r2c);
    pm_r2c(pm, canvas, delta_k);
    LEAVE(r2c);

    /* calculate the forces save them to p->acc */
    apply_dealiasing_transfer(gravity, pm, delta_k, delta_k);


    /*emit before force event*/
    CLOCK(beforeforce);
    FastPMForceEvent eventb[1];

    FastPMPowerSpectrum pscdm[1];
    FastPMPowerSpectrum psnu[1];
    double fnu = 0.006;
    // Or only current time step read?
    char filecdm[1024];
    char filenu[1024];

    sprintf(filecdm, "/home/energy/grad_fPM/cleanPM/fastpm/tests/cdm%1.3f.dat", trans->a.f);
    sprintf(filenu, "/home/energy/grad_fPM/cleanPM/fastpm/tests/nu%1.3f.dat", trans->a.f);
//    printf("cdm filename=%s\n", file);
    printf("steps = %1.3f\n",trans->a.f);
    char *contentcdm = fastpm_file_get_content(filecdm);//casting?!
    char *contentnu = fastpm_file_get_content(filenu);//casting?!

//    if(content!=0) printf("Content%s!\n",content);
    printf("Hi hi 1\n");
    fastpm_powerspectrum_init_from_string(pscdm, contentcdm);
    fastpm_powerspectrum_init_from_string(psnu, contentnu);
    printf("Hi hi ih\n");
    free(contentcdm);
    free(contentnu);
    //    printf("Hi hi! %g",fastpm_powerspectrum_eval(ps,0.001));
    eventb->step += 1;
    eventb->a_f = trans->a.f;
    PMKIter kiter;
    double k0 = 2 * M_PI / pm_boxsize(pm)[0];
    for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
        int d;
        ptrdiff_t kk = 0.;
        ptrdiff_t ind = kiter.ind;
        for(d = 0; d < 3; d++) {
            double ik = kiter.iabs[d];
            if(ik > pm->Nmesh[d] / 2) ik -= pm->Nmesh[d];
            kk += ik * ik;
        }
        double k = sqrt(kk) * k0;
//        printf("k is %f\n",k);
//        printf("transfer is %f\n",ps->k[1]);
//        if(fastpm_powerspectrum_eval(ps, k)) printf("evaluated %g\n",fabs(fastpm_powerspectrum_eval(ps, k)));
        double Transfer = ((1-fnu)*fabs(fastpm_powerspectrum_eval(pscdm, k))
                            +fnu*fabs(fastpm_powerspectrum_eval(pscdm, k)))
                            / fabs(fastpm_powerspectrum_eval(pscdm, k));
//        delta_k[ind+0] = delta_k[ind+0]*Transfer;
//        delta_k[ind+1] = delta_k[ind+1]*Transfer;
        eventb->delta_k = delta_k;
    }

    fastpm_solver_emit_event(fastpm, FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_BEFORE, (FastPMEvent*) eventb);
    LEAVE(beforeforce);

    int d;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z, PACK_POTENTIAL};
    for(d = 0; d < (gravity->ComputePotential?4:3); d ++) {
        CLOCK(transfer);
        apply_kernel_transfer(gravity, pm, delta_k, canvas, d);
        LEAVE(transfer);

        CLOCK(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        CLOCK(readout);
        fastpm_readout_local(reader, canvas, p, p->np + pgd->nghosts, NULL, ACC[d]);
        LEAVE(readout);

        CLOCK(reduce);
        pm_ghosts_reduce(pgd, ACC[d]);
        LEAVE(reduce);
    }

    CLOCK(afterforce);
    FastPMForceEvent event[1];
    event->delta_k = delta_k;
    event->a_f = trans->a.f;

    fastpm_solver_emit_event(fastpm, FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER, (FastPMEvent*) event);
    LEAVE(afterforce);


    pm_free(pm, canvas);
    pm_free(pm, delta_k);
    pm_ghosts_free(pgd);
}

void
fastpm_gravity_calculate_gradient(FastPMGravity * gravity,
    PM * pm,
    FastPMStore * grad_acc,
    FastPMStore * p,
    FastPMStore * grad_pos
)
{
    FastPMFloat * grad_delta_k = pm_alloc(pm);
    FastPMFloat * grad_canvas = pm_alloc(pm);
    FastPMFloat * canvas = pm_alloc(pm);
    FastPMFloat * delta_k = pm_alloc(pm);

    FastPMStore grad_pos_1[1];

    fastpm_store_init(grad_pos_1, grad_pos->np_upper, grad_pos->attributes);

    grad_pos_1->np = grad_pos->np;

    FastPMPainter reader[1];
    FastPMPainter painter[1];

    fastpm_painter_init(reader, pm, gravity->PainterType, gravity->PainterSupport);
    fastpm_painter_init(painter, pm, gravity->PainterType, gravity->PainterSupport);

    /* watch out: boost the density since mesh is finer than grid */
    long long np = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, pm_comm(pm));

    double density_factor = pm->Norm / np;

    CLOCK(paint);
    fastpm_paint(painter, canvas, p, NULL, 0);
    fastpm_apply_multiply_transfer(pm, canvas, canvas, density_factor);
    LEAVE(paint);

    CLOCK(r2c);
    pm_r2c(pm, canvas, delta_k);
    LEAVE(r2c);

    /* calculate the forces save them to p->acc */
    apply_dealiasing_transfer(gravity, pm, delta_k, delta_k);

    int d;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};

    memset(grad_delta_k, 0, sizeof(grad_delta_k[0]) * pm_allocsize(pm));
    memset(grad_canvas, 0, sizeof(grad_canvas[0]) * pm_allocsize(pm));
    memset(&grad_pos->x[0][0], 0, sizeof(grad_pos->x[0]) * grad_pos->np_upper);

    for(d = 0; d < 3; d ++) {
        FastPMFloat * grad_delta_k_1 = pm_alloc(pm);
        FastPMFloat * grad_canvas_1 = pm_alloc(pm);

        /* gradient of read out */
        CLOCK(transfer);
        apply_kernel_transfer(gravity, pm, delta_k, canvas, d);
        pm_assign(pm, delta_k, canvas);
        LEAVE(transfer);

        CLOCK(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        CLOCK(readout);
        fastpm_readout_gradient(reader, grad_acc,
            canvas, p, NULL, ACC[d], grad_canvas_1, grad_pos_1);
        LEAVE(readout);

        pm_c2r_gradient(pm, grad_canvas_1, grad_delta_k_1);
        apply_kernel_transfer(gravity, pm, grad_delta_k_1, grad_delta_k_1, d);

        ptrdiff_t i;
        for(i = 0; i < pm_allocsize(pm); i ++) {
            grad_delta_k[i] += grad_delta_k_1[i];
        }

        for(i = 0; i < grad_pos->np; i ++) {
            int d;
            for(d = 0; d < 3; d ++) {
                grad_pos->x[i][d] += grad_pos_1->x[i][d];
            }
        }

        pm_free(pm, grad_canvas_1);
        pm_free(pm, grad_delta_k_1);
    }

    pm_r2c_gradient(pm, grad_delta_k, grad_canvas);

    apply_dealiasing_transfer(gravity, pm, grad_delta_k, grad_delta_k);

    fastpm_paint_gradient(painter, grad_canvas, p, NULL, 0, grad_pos_1, NULL);

    ptrdiff_t i;
    for(i = 0; i < grad_pos->np; i ++) {
        int d;
        for(d = 0; d < 3; d ++) {
            grad_pos->x[i][d] += density_factor * grad_pos_1->x[i][d];
        }
    }

    fastpm_store_destroy(grad_pos_1);
    pm_free(pm, delta_k);
    pm_free(pm, canvas);
    pm_free(pm, grad_canvas);
    pm_free(pm, grad_delta_k);
}
