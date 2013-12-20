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
#include "parseutils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void parseutils_buf2str(char *out, const char *buf, size_t len)
{
    memcpy(out, buf, len);
    out[len] = '\0';
}
int parseutils_buf2int(int *out, const char *buf, size_t len)
{
    char tmp[20];
    if (len>=sizeof(tmp)) return -1;
    parseutils_buf2str(tmp, buf, len);
    *out = strtol(tmp, NULL, 10);
    return 0;
}
int parseutils_buf2double(double *out, const char *buf, size_t len)
{
    char tmp[30];
    if (len>=sizeof(tmp)) return -1;
    parseutils_buf2str(tmp, buf, len);
    *out = strtod(tmp, NULL);
    return 0;
}
/* decimal point assumed */
int parseutils_buf2doubleDPA(double *out, const char *buf, size_t len)
{
    char tmp[30], newT[50] = {0};
    int i, j;
    if (len>=sizeof(tmp)) return -1;
    parseutils_buf2str(tmp, buf, len);
    i = j = 0;
    while (tmp[i]==' ') ++i;
    if (tmp[i]=='-') {
        strcat(newT, "-0.");
        j = 3;
        ++i;
    } else {
        strcat(newT, "0.");
        j = 2;
    }
    while (tmp[i]>='0' && tmp[i]<='9') {
        newT[j++] = tmp[i++];
        newT[j] = '\0';
    }
    if (tmp[i]=='-' || tmp[i]=='+') {
        strcat(newT, "e");
        strcat(newT, tmp+i);
    }
    *out = strtod(newT, NULL);
    return 0;
}

