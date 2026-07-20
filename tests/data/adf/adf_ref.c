/*
 * adf_ref.c - native parity reference for the PQAnalysis ADF analysis.
 *
 * This is a standalone, double-precision port of the legacy thh_tools
 * ADF algorithm that survives as commented code in
 * RDF/source/process.c (lines ~210-255). For every center atom i, every
 * ligand-1 atom j (j != i) and every ligand-2 atom k (k != i, k != j)
 * it computes the orthorhombic minimum image vectors v_ij and v_ik with
 * the legacy convention
 *
 *     x2   = CRD[j] + BOX * round((CRD[i] - CRD[j]) / BOX);
 *     dx12 = CRD[i] - x2;                 // = image vector j -> i
 *
 * the i-j and i-k distances r12, r13, applies the optional radial gates
 * (r < r_min or r >= r_max skips) and bins the j-i-k angle
 *
 *     angle = 180/PI * acos((v_ij . v_ik) / r12 / r13);
 *     bin   = (int) floor(angle / delta_angle);
 *
 * Ordered ligand pairs are counted, exactly as the legacy tool did.
 *
 * Two documented deviations from the raw legacy remnant make the
 * reference robust (and match the PQAnalysis port): the acos argument is
 * clamped to [-1, 1] before acos (the legacy code did not, a latent NaN
 * bug) and an exact 180 degree triplet, which lands in bin n_bins, is
 * discarded instead of written out of bounds (the legacy code would have
 * written past the end of adf_hist).
 *
 * Usage:
 *   adf_ref TRAJ.xyz DELTA_ANGLE CENTER LIG1 LIG2 \
 *           [R_MIN_1 R_MAX_1 R_MIN_2 R_MAX_2]
 *
 * CENTER/LIG1/LIG2 are element names selected from the xyz atom column.
 * A gate is disabled by passing a negative bound. The histogram is
 * printed as "bin_center count" lines (n_bins = floor(180/DELTA_ANGLE)).
 *
 * Build: cc -O2 -o adf_ref adf_ref.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_ATOMS 100000

int main(int argc, char **argv)
{
    if (argc != 6 && argc != 10)
    {
        fprintf(stderr,
                "usage: %s TRAJ DELTA CENTER LIG1 LIG2 "
                "[RMIN1 RMAX1 RMIN2 RMAX2]\n",
                argv[0]);
        return 1;
    }

    const char *traj_file = argv[1];
    double delta_angle = atof(argv[2]);
    const char *center_name = argv[3];
    const char *lig1_name = argv[4];
    const char *lig2_name = argv[5];

    int gate_1 = 0, gate_2 = 0;
    double r_min_1 = 0.0, r_max_1 = 0.0, r_min_2 = 0.0, r_max_2 = 0.0;

    if (argc == 10)
    {
        r_min_1 = atof(argv[6]);
        r_max_1 = atof(argv[7]);
        r_min_2 = atof(argv[8]);
        r_max_2 = atof(argv[9]);
        /* a negative bound disables that side; a gate is active if
         * either of its bounds is non-negative */
        gate_1 = (r_min_1 >= 0.0) || (r_max_1 >= 0.0);
        gate_2 = (r_min_2 >= 0.0) || (r_max_2 >= 0.0);
        if (r_min_1 < 0.0) r_min_1 = 0.0;
        if (r_max_1 < 0.0) r_max_1 = INFINITY;
        if (r_min_2 < 0.0) r_min_2 = 0.0;
        if (r_max_2 < 0.0) r_max_2 = INFINITY;
    }

    int n_bins = (int) floor(180.0 / delta_angle);
    if (n_bins < 1)
    {
        fprintf(stderr, "delta_angle too large\n");
        return 1;
    }

    long *hist = (long *) calloc((size_t) n_bins, sizeof(long));

    FILE *fp = fopen(traj_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "cannot open %s\n", traj_file);
        return 1;
    }

    static double crd[MAX_ATOMS][3];
    static char name[MAX_ATOMS][16];
    static int is_center[MAX_ATOMS];
    static int is_lig1[MAX_ATOMS];
    static int is_lig2[MAX_ATOMS];

    char line[4096];
    const double rad2deg = 180.0 / M_PI;

    while (fgets(line, sizeof(line), fp) != NULL)
    {
        int n_atoms = 0;
        double bx = 0.0, by = 0.0, bz = 0.0;

        /* header line: n_atoms bx by bz (orthorhombic) */
        if (sscanf(line, "%d %lf %lf %lf", &n_atoms, &bx, &by, &bz) < 4)
        {
            /* allow a header without a box (single value) -> skip */
            if (sscanf(line, "%d", &n_atoms) != 1)
                continue;
            bx = by = bz = INFINITY;
        }

        /* comment line */
        if (fgets(line, sizeof(line), fp) == NULL)
            break;

        for (int a = 0; a < n_atoms; a++)
        {
            if (fgets(line, sizeof(line), fp) == NULL)
            {
                fclose(fp);
                free(hist);
                fprintf(stderr, "unexpected end of file\n");
                return 1;
            }

            sscanf(line, "%15s %lf %lf %lf",
                   name[a], &crd[a][0], &crd[a][1], &crd[a][2]);

            is_center[a] = (strcmp(name[a], center_name) == 0);
            is_lig1[a] = (strcmp(name[a], lig1_name) == 0);
            is_lig2[a] = (strcmp(name[a], lig2_name) == 0);
        }

        double box[3] = {bx, by, bz};

        for (int i = 0; i < n_atoms; i++)
        {
            if (!is_center[i])
                continue;

            for (int j = 0; j < n_atoms; j++)
            {
                if (!is_lig1[j] || j == i)
                    continue;

                double v12[3], v13[3];
                double r12sq = 0.0;

                for (int d = 0; d < 3; d++)
                {
                    double x2 = crd[j][d]
                        + box[d] * round((crd[i][d] - crd[j][d]) / box[d]);
                    v12[d] = crd[i][d] - x2;
                    r12sq += v12[d] * v12[d];
                }

                double r12 = sqrt(r12sq);

                if (gate_1 && (r12 < r_min_1 || r12 >= r_max_1))
                    continue;

                for (int k = 0; k < n_atoms; k++)
                {
                    if (!is_lig2[k] || k == i || k == j)
                        continue;

                    double r13sq = 0.0;

                    for (int d = 0; d < 3; d++)
                    {
                        double x3 = crd[k][d]
                            + box[d]
                                  * round((crd[i][d] - crd[k][d]) / box[d]);
                        v13[d] = crd[i][d] - x3;
                        r13sq += v13[d] * v13[d];
                    }

                    double r13 = sqrt(r13sq);

                    if (gate_2 && (r13 < r_min_2 || r13 >= r_max_2))
                        continue;

                    double dot = v12[0] * v13[0]
                        + v12[1] * v13[1]
                        + v12[2] * v13[2];

                    double cosine = dot / r12 / r13;

                    if (cosine > 1.0)
                        cosine = 1.0;
                    else if (cosine < -1.0)
                        cosine = -1.0;

                    double angle = rad2deg * acos(cosine);

                    int bin = (int) floor(angle / delta_angle);

                    if (bin >= 0 && bin < n_bins)
                        hist[bin] += 1;
                }
            }
        }
    }

    fclose(fp);

    for (int b = 0; b < n_bins; b++)
        printf("%.4f %ld\n", ((double) b + 0.5) * delta_angle, hist[b]);

    free(hist);

    return 0;
}
