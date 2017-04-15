FASTPM_BEGIN_DECLS

void
fastpm_solver_gravity_calculate(FastPMSolver * fastpm, FastPMGravity * gravity, PM * pm, FastPMStore * p, FastPMTransition * trans);

void
fastpm_gravity_calculate_gradient(FastPMGravity * gravity,
        PM * pm,
        FastPMStore * grad_acc,
        FastPMStore * p,
        FastPMStore * grad_pos
        );

FASTPM_END_DECLS
