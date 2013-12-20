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
#include "tle.h"
#include "parseutils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int computeChk(const char *buf)
{
    int i, chk=0;
    for (i=0; i<68; ++i) {
        if (buf[i]>='0' && buf[i]<='9')
            chk += buf[i]-'0';
        if (buf[i]=='-')
            chk += 1;
    }
    return chk%10;
}

static int header_parse(tle_st *tle, const char *buf)
{
    size_t len = strlen(buf);
    int i;
    while (buf[len-1]=='\n' || buf[len-1]=='\r') --len;
    if (len==0 || len>=sizeof(tle->name))
        return -1;
    strcpy(tle->name, buf);
    for (i=strlen(tle->name)-1; i>=0; --i)
        if (tle->name[i]=='\r' || tle->name[i]=='\n' || tle->name[i]==' ') {
            tle->name[i]='\0';
            continue;
        } else
            break;
    if (i<=0) return -1;
    return 0;
}

static int line1_parse(tle_st *tle, const char *buf)
{
    if (strlen(buf)<69) return -1;
    if (buf[0]!='1') return -1;
    if (parseutils_buf2int(&tle->satNo, buf+2, 5)!=0) return -1;
    tle->class = buf[7];
    if (parseutils_buf2int(&tle->designatorYY, buf+9, 2)!=0) return -1;
    if (parseutils_buf2int(&tle->designatorLaunchNo, buf+11, 3)!=0) return -1;
    parseutils_buf2str(tle->designatorPiece, buf+14, 3);
    if (parseutils_buf2int(&tle->epochYY, buf+18, 2)!=0) return -1;
    if (parseutils_buf2double(&tle->epoch, buf+20, 12)!=0) return -1;
    if (parseutils_buf2double(&tle->dtmm, buf+33, 10)!=0) return -1;
    if (parseutils_buf2doubleDPA(&tle->dt2mm, buf+44, 8)!=0) return -1;
    if (parseutils_buf2doubleDPA(&tle->bstar, buf+53, 8)!=0) return -1;
    /*if (buf[62]!='0') return -1;*/
    if (parseutils_buf2int(&tle->elSet, buf+64, 4)!=0) return -1;
    if (computeChk(buf)!=buf[68]-'0') {
        fprintf(stderr, "checksum error: got %d expected %d\n",
                buf[68]-'0', computeChk(buf));
        return -1;
    }
    return 0;
}

static int line2_parse(tle_st *tle, const char *buf)
{
    int sN;
    if (strlen(buf)<69) return -1;
    if (buf[0]!='2') return -1;
    if (parseutils_buf2int(&sN, buf+2, 5)!=0 || sN!=tle->satNo) return -1;
    if (parseutils_buf2double(&tle->inclinationDeg, buf+8, 8)!=0) return -1;
    if (parseutils_buf2double(&tle->rightAscDeg, buf+17, 8)!=0) return -1;
    if (parseutils_buf2doubleDPA(&tle->ecc, buf+26, 7)!=0) return -1;
    if (parseutils_buf2double(&tle->argOfPerigeeDeg, buf+34, 8)!=0) return -1;
    if (parseutils_buf2double(&tle->meanAnomalyDeg, buf+43, 8)!=0) return -1;
    if (parseutils_buf2double(&tle->meanMotionRpD, buf+52, 11)!=0) return -1;
    if (parseutils_buf2int(&tle->revs, buf+63, 5)!=0) return -1;
    if (computeChk(buf)!=buf[68]-'0') {
        fprintf(stderr, "checksum error: got %d expected %d\n",
                buf[68]-'0', computeChk(buf));
        return -1;
    }
    return 0;
}

void tle_destroy(tle_st **tle)
{
    tle_st *p;
    if (tle==NULL || *tle==NULL) return;
    p = (*tle)->next;
    free(*tle);
    *tle = NULL;
    tle_destroy(&p);
}

tle_st *tle_read(FILE *in, int readAll, int compatMode)
{
    char buf[150];
    tle_st *p, *head;
    int recpos, line;
    head = NULL;
    recpos = 0;
    line = 1;
    while (fgets(buf, sizeof(buf), in)!=NULL) {
        if (buf[strlen(buf)-1]!='\n' && buf[strlen(buf)-1]!='\r') goto badexit;
        if (compatMode && buf[0]=='#') {
            ++line;
            continue;
        }
        if (recpos==0) {
            if ((p = malloc(sizeof(*p)))==NULL)
                goto badexit;
            p->next = head;
            head = p;
            if (compatMode) {
                ++recpos;
                goto line1parse;
            }
            if (header_parse(p, buf)!=0) goto badexit;
            ++recpos;
        } else if (recpos==1) {
            line1parse:
            if (line1_parse(p, buf)!=0) goto badexit;
            ++recpos;
        } else if (recpos==2) {
            if (line2_parse(p, buf)!=0) goto badexit;
            recpos = 0;
            if (!readAll)
                break;
        }
        ++line;
    }
    return head;

badexit:
    fprintf(stderr, "error parsing file at line %d\n", line);
    if (p!=NULL)
        fprintf(stderr, "\tlast sat #: %d\n", p->satNo);
    tle_destroy(&head);
    return NULL;
}

/*int main(void)
{
    tle_st *p;
    p = tle_read("visual.txt", 1, 0);
    printf((p==NULL)?"nope\n":"ok\n");
    tle_destroy(&p);
    return 0;
}
*/
