#include <stdio.h>
#include <string.h>
#include <alloca.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    int ThisTask;
    int NTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);

    FastPMColumn x[1];
    int64_t index[128];
    fastpm_column_init_float3(x, 256);

    size_t oldsize = (1 + ThisTask) * 10;
    size_t newsize = (NTask - ThisTask) * 10;
    fastpm_column_resize(x, oldsize);
    int i;

    for(i = 0; i < oldsize; i ++) {
        double pos[3] = {ThisTask * 10 + i, ThisTask * 10 + i, i};
        fastpm_column_set_double(x, i, pos);
        index[i] = i;
    }

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("oldsize = %td\n", oldsize);
    for(i = 0; i < oldsize; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("oldpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_parallel_permute(x, index, newsize, comm);

    for(i = 0; i < ThisTask; i ++) {
        MPI_Barrier(comm);
    }
    printf("newsize = %td\n", newsize);
    for(i = 0; i < newsize; i ++) {
        double pos[3];
        fastpm_column_get_double(x, i, pos);
        printf("newpos[%d] =%g %g %g\n", i, pos[0], pos[1], pos[2]);
    }
    for(i = ThisTask; i < NTask; i ++) {
        MPI_Barrier(comm);
    }

    fastpm_column_destroy(x);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
