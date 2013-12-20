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
#include "vsop87file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parseutils.h"

/* file extension */
static const char *s_ext[] = {
    "",
    "mer", "ven", "ear", "mar", "jup", "sat", "ura", "nep", "emb", "sun"
};

static int header_verify(vsop87_hdr_st *hdr)
{
    if (hdr->versionCode<0 || hdr->versionCode>5) return -1;
    if (hdr->coordIndex<1 || hdr->coordIndex>6) return -1;
    if (hdr->nTerms<=0) return -1;
    return 0;
}
static int header_parse(vsop87_hdr_st *hdr, const char *buf)
{
    hdr->versionCode = buf[17]-'0';
    parseutils_buf2str(hdr->name, buf+22, 7);
    hdr->coordIndex = buf[41]-'0';
    hdr->alpha = buf[59]-'0';
    if (parseutils_buf2int(&hdr->nTerms, buf+60, 7)!=0) return -1;
    return header_verify(hdr);
}

static int record_verify(vsop87_rcd_st *rcd)
{
    if (rcd->versionCode<0 || rcd->versionCode>5) return -1;
    if (rcd->bodyCode<1 || rcd->bodyCode>9) return -1;
    if (rcd->coordIndex<1 || rcd->coordIndex>6) return -1;
    return 0;
}
static int record_parse(vsop87_rcd_st *rcd, const char *buf)
{
    int i;
    rcd->versionCode = buf[1]-'0';
    rcd->bodyCode = buf[2]-'0';
    rcd->coordIndex = buf[3]-'0';
    rcd->alpha = buf[4]-'0';
    if (parseutils_buf2int(&rcd->n, buf+5, 5)!=0) return -1;
    for (i=0; i<12; ++i)
        if (parseutils_buf2int(&rcd->a[i], buf+10+i*3, 3)!=0) return -1;
    if (parseutils_buf2double(&rcd->S, buf+46, 15)!=0) return -1;
    if (parseutils_buf2double(&rcd->K, buf+61, 18)!=0) return -1;
    if (parseutils_buf2double(&rcd->A, buf+79, 18)!=0) return -1;
    if (parseutils_buf2double(&rcd->B, buf+97, 14)!=0) return -1;
    if (parseutils_buf2double(&rcd->C, buf+111, 20)!=0) return -1;
    return record_verify(rcd);
}

void vsop87_filepart_destroy(vsop87_filepart_st **fp)
{
    vsop87_filepart_st *p;
    if (fp==NULL || *fp==NULL) return;
    p = (*fp)->next;
    free((*fp)->rcds);
    free(*fp);
    *fp = NULL;
    vsop87_filepart_destroy(&p);
}

vsop87_filepart_st *vsop87_filepart_read(const char *datadir,
                                         vsop87_body_et body, char code)
{
    char *fname;
    char buf[150];
    vsop87_filepart_st *p, *head;
    int recpos, line;
    FILE *in;
    int addSlash = (datadir[strlen(datadir)-1]!='/');
    if ((fname = malloc(strlen(datadir)+addSlash+12))==NULL)
        return NULL;
    sprintf(fname, "%s", datadir);
    if (addSlash) strcat(fname, "/");
    if (code!=' ') {
        if (body==BODY_EARTHMOON && code=='E')
            sprintf(fname, "%sVSOP87%c.sun", fname, code);
        else
            sprintf(fname, "%sVSOP87%c.%s", fname, code, s_ext[body]);
    }
    else sprintf(fname, "%sVSOP87.%s", fname, s_ext[body]);
    if ((in = fopen(fname, "rb"))==NULL) {
        perror("fopen");
        free(fname);
        return NULL;
    }
    head = NULL;
    recpos = 0;
    line = 1;
    while (fgets(buf, sizeof(buf), in)!=NULL) {
        if (buf[strlen(buf)-1]!='\n' && buf[strlen(buf)-1]!='\r') goto badexit;
        if (recpos==0) {
            if ((p = malloc(sizeof(*p)))==NULL)
                goto badexit;
            p->rcds = NULL;
            p->next = head;
            head = p;
            if (header_parse(&p->hdr, buf)!=0) goto badexit;
            if ((p->rcds = calloc(p->hdr.nTerms, sizeof(*p->rcds)))==NULL)
                goto badexit;
            ++recpos;
        } else {
            if (record_parse(&p->rcds[recpos-1], buf)!=0) goto badexit;
            ++recpos;
            if (recpos>p->hdr.nTerms)
                recpos = 0;
        }
        ++line;
    }
    fclose(in);
    free(fname);
    return head;

badexit:
    printf("error parsing file %s at line %d\n", fname, line);
    fclose(in);
    free(fname);
    vsop87_filepart_destroy(&head);
    return NULL;
}

/*int main(void)
{
    vsop87_filepart_st *ret;
    ret = vsop87_filepart_read(BODY_EARTH);
    printf((ret==NULL)?"nope\n":"ok\n");
    vsop87_filepart_destroy(&ret);
    return 0;
}*/

