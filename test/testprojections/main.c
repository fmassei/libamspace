#include <libamspace.h>
#include <stdio.h>
#include <stdlib.h>
#include <projections.h>
#include <angles.h>

/* general vector checker */
#define ISNEAR(_X, _Y) (FABS((_X)-(_Y))<1e-5)
static int matchRes(DT *v1, DT *v2, int n)
{
    int i;
    for (i=0; i<n; ++i)
        if (!ISNEAR(v1[i], v2[i]))
            return 0;
    return 1;
}
#undef ISNEAR
int main(void)
{
    DT r[][2] = {
        {  60, -120 },
        {  30,  -88 },
        {  17,  -22 },
        {  30,   77 },
        {  60,  121 },
        {  40,  149 },
        { -38,  133 },
        { -44,   32 },
        { -77,  142 },
        { -65, -142 },
        { -10,  -20 },
        { -55,  -79 },
        {   0,    3 },
        {   0,    0 },
    };
    DT res[][3] = {
        { 6655207.07100948,    332705.152498051,    11 },
        { 3319206.45938262,    403549.847394392,    16 },
        {  1879826.8019149,    393552.252927449,    27 },
        { 3320470.02816103,    692915.106411067,    43 },
        { 6653098.17877781,    388455.954548937,    51 },
        { 4429673.81714971,    670725.499572411,    55 },
        { 5792296.80109494,    324396.624426642,    53 },
        { 5127640.84111201,    419825.73997373,     36 },
        {  1453013.1989751,    525110.200746697,    54 },
        { 2791172.24442381,    452844.877550694,     7 },
        { 8894421.31703289,    609600.772154029,    27 },
        { 3903378.48696447,    627928.196137513,    17 },
        {                0,              500000,    31 },
        {                0,    166021.549696835,    31 }
    };
    DT out[4];
    gravconst_st g;
    int i;
    getgravconst(WGS84, &g);
    for (i=0; i<14; ++i) {
        amsp_projections_latlon2utm(out, r[i], &g);
        if (!matchRes(out, res[i], 3)) {
            printf("%d failed (%f %f %f)\n", i, out[0], out[1], out[2]);
            return EXIT_FAILURE;
        }
    }
    return 0;
}

