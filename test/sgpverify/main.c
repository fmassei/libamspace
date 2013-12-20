/*
    This file is part of libamspace.

    libamspace is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libamspace is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libamspace.h>

int main(const int argc, const char *argv[])
{
    FILE *in1, *in2;
    char buf1[400], buf2[400];
    double sum_ep, sum_ev;
    int sum_n;
    int line = 0, summarize=0;

    if (argc<2) {
        fprintf(stderr, "Usage: sgpverify <file1|-> <file2>\n");
        return EXIT_FAILURE;
    }
    if (!strcmp(argv[1], "-")) {
        in1 = stdin;
    } else {
        if ((in1 = fopen(argv[1], "rb"))==NULL) {
            perror("fopen");
            return EXIT_FAILURE;
        }
    }
    if ((in2 = fopen(argv[2], "rb"))==NULL) {
        perror("fopen");
        return EXIT_FAILURE;
    }
    if (argc>3)
        summarize = 1;
    while ( fgets(buf1, sizeof(buf1), in1)!=NULL &&
            fgets(buf2, sizeof(buf2), in2)!=NULL) {
        double d1, d2, val1[6], val2[6];
        double ep, ev;
        int j;
        ++line;
        if (    (buf1[strlen(buf1)-1]!='\n' && buf1[strlen(buf1)-1]!='\r') ||
                (buf2[strlen(buf2)-1]!='\n' && buf2[strlen(buf2)-1]!='\r')  ) {
            fprintf(stderr, "unexpected long line %d!\n", line);
            return EXIT_FAILURE;
        }
        /* check for sat header */
        if (strstr(buf1, "xx")!=NULL) {
            if (strstr(buf2, "xx")==NULL) {
                fprintf(stderr, "test data disalignment on line %d!\n", line);
                return EXIT_FAILURE;
            }
            sscanf(buf1, "%lf", &d1);
            sscanf(buf2, "%lf", &d2);
            if (d1!=d2) {
                fprintf(stderr, "test data disalignment on line %d!\n", line);
                return EXIT_FAILURE;
            }
            if (summarize) {
                if (line!=1)
                    printf("%lf %lf\n", sum_ep/sum_n, sum_ev/sum_n);
                sum_ep = sum_ev = sum_n = 0.;
            }
            printf("\n%lf xx\n", d1);
            continue;
        }
        sscanf(buf1, "%lf %lf %lf %lf %lf %lf %lf",
            &d1, &val1[0], &val1[1], &val1[2], &val1[3], &val1[4], &val1[5]);
        sscanf(buf2, "%lf %lf %lf %lf %lf %lf %lf",
            &d2, &val2[0], &val2[1], &val2[2], &val2[3], &val2[4], &val2[5]);
        if (d1!=d2) {
            fprintf(stderr, "test data disalignment on line %d!\n", line);
            return EXIT_FAILURE;
        }
        ep = SQRT(P2(val1[0]-val2[0])+
                  P2(val1[1]-val2[1])+
                  P2(val1[2]-val2[2]));
        ev = SQRT(P2(val1[3]-val2[3])+
                  P2(val1[4]-val2[4])+
                  P2(val1[5]-val2[5]));
        if (summarize) {
            sum_ep += ep;
            sum_ev += ev;
            ++sum_n;
        } else {
            for (j=0; j<3; ++j)
                printf("%lf ", val2[j]-val1[j]);
            printf("[%lf] ", ep);
            for (j=3; j<6; ++j)
                printf("%lf ", val2[j]-val1[j]);
            printf("[%lf] ", ev);
            printf("\n");
        }
    }
    if (summarize)
        printf("%lf %lf\n", sum_ep/sum_n, sum_ev/sum_n);

    fclose(in2);
    fclose(in1);
    return 0;
}

