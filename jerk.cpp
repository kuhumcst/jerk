/* JERK
Analyses the movement of two points in x-y plane, in casu nose tips data from
OpenPoseDemo.exe, and computes velocity, acceleration and jerk of the points.

The program expects input with data for two persons in a dialogue situation.

The format of the input is as follows:

FrameNr <tab> <person A data> <tab> <person B data> <newline>
  ...

where data for each person is

<x-position> <tab> <y-position> <tab> <probability>

The probability value is currently ignored.

The command line arguments to the program are
  input file
  output file
  "name" of the first person
  "name" of the second person
  number of observations used for computing velocity
  number of observations used for computing acceleration
  number of observations used for computing jerk
  location of the observations with respect to the created annotation:
       'past' or 'middle'

For example

  jerk F2_M4_SPLIT_FINAL.tab F2_M4_SPLIT_FINAL-7-14-21-past.vaj.tab F2 M4 7 14 21 past


Output format is
<frame number x 4> <tab> <person A data> <tab> <person B data> <new line>
  ...

<person data> =
<x> <tab> <y> <tab> <probability> <tab> <velocity> <tab> <acceleration> <tab> <jerk> <new line>

where
  <velocity> = <strength> <tab> <angle (1..12)> <tab> <x component> <tab> <y component>
   etc.

<strength> and <angle> are just another (polar) representation of the components.
1..12 are 'clock' directions. '12' is up, '3' to the right, etc.

Author: Bart Jongejan
2020.01.23
*/

/*  */

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define M_PI           3.14159265358979323846  /* pi */

enum
    {
    evelo, eacce, ejerk
    };

struct person
    {
    double x;
    double y;
    double p;
    double vx;
    double ax;
    double jx;
    double vy;
    double ay;
    double jy;
    };

int persons;

struct line
    {
    int frame;
    struct person* P;
    line() :frame(-1), P(0)
        {
        P = new person[persons];
        }
    };

double velocity(double St, double St2, double Sh, double Sth, double period)
    {
    /*
    first degree: a:position, b:velocity
    var a
    solution (a.period^-1*(Sh+-1*St*b))
    var b
      solution
      (b.(-1*St^2+St2*period)^-1*(-1*Sh*St+Sth*period))

    For vertical, replace h th by v tv
    */
    return (Sth * period - Sh * St) / (St2 * period - St * St);
    }

double acceleration(double St, double St2, double St3, double St4, double Sh, double Sth, double St2h, double period)
    {

    /*

    acceleration c

    (horizontal)

    accumulate (S) t t2 t3 t4 h th t2h
    (vertical                 v tv t2v )

    solution
    ( c
    .     ( -1*St2^3
        + 2*St*St2*St3
        + St2*St4*period
        + -1*St^2*St4
        + -1*St3^2*period
        )
      ^ -1
    * ( -1*Sh*St2^2
      + Sh*St*St3
      + St*St2*Sth
      + St2*St2h*period
      + -1*St3*Sth*period
      + -1*St^2*St2h
      )
    )

    For vertical, replace h th t2h by v tv t2v
    */
    return
        (Sh * (St * St3 - St2 * St2)
            + Sth * (St * St2 - St3 * period)
            + St2h * (St2 * period - St * St)
            )
        /
        (St2 * (2 * St * St3 - St2 * St2 + St4 * period)
            - St * St * St4
            - St3 * St3 * period
            );

    }

double jerk(double St2, double St2h, double St3, double St3h, double St4, double St5, double St6, double Sth, double period)
    {
    double var1 = St3 * St3;
    double var2 = St4 * St4;
    double var5 = St3 * St5;
    double var6 = var1 - St2 * St4;

    return -(St2h * (St2 * St5 * period - St3 * (St2 * St2 + St4 * period)) + St3h * (St2 * St2 * St2 + period * var6) + Sth * (St2 * var6 + period * (var2 - var5)))
        / (var1 * var1 - St2 * (St5 * St5 * period - St2 * (var2 + 2 * var5 - St2 * St6) + St4 * (3 * var1 - St6 * period)) - period * (St6 * var1 + St4 * (var2 + -2.0 * var5)));
    }


void doPerson(int Person, int what, const char* Period, const char* margin, line* Lines)
    {
    size_t seqsiz = strtoul(Period, 0, 10);
    double period = (double)seqsiz;
    long t0;
    long* Ts = (long*)calloc(seqsiz, sizeof(long));
    double* ts = (double*)calloc(seqsiz, sizeof(double));
    double* hs = (double*)calloc(seqsiz, sizeof(double));
    double* vs = (double*)calloc(seqsiz, sizeof(double));
    size_t index = 0;
    struct line* current = Lines;

    long Offset;
    if (!strcmp(margin, "middle"))
        Offset = (seqsiz - 1) >> 1;
    else
        Offset = 0;
    for (current = Lines; current->frame != 0; ++current)
        {
        struct person* pers = current->P + Person;
        struct line* assignTo = current - Offset;
        struct person* assignToPers = assignTo->P + Person;
        if (assignTo < Lines)
            assignTo = 0;
        else
            switch (what)
                {
                case evelo:
                    assignToPers->vx = 0;
                    assignToPers->vy = 0;
                    break;
                case eacce:
                    assignToPers->ax = 0;
                    assignToPers->ay = 0;
                    break;
                case ejerk:
                    assignToPers->jx = 0;
                    assignToPers->jy = 0;
                    break;
                }

        if (pers->x != 0)
            {
            int ind = index % seqsiz;

            if (index == 0)
                // Start accumulation of datapoints
                {
                t0 = current->frame;
                Ts[ind] = 0;
                }
            else
                {
                Ts[ind] = current->frame - t0;
                }
            hs[ind] = pers->x;
            vs[ind] = pers->y;
            ++index;
            if (index >= seqsiz)
                {
                size_t i;
                double averageTime = 0;

                double St = 0.0, St2 = 0.0, St3 = 0.0, St4 = 0.0, St5 = 0.0, St6 = 0.0,
                    Sh = 0.0, Sth = 0.0, St2h = 0.0, St3h = 0.0,
                    Sv = 0.0, Stv = 0.0, St2v = 0.0, St3v = 0.0;
                double averageh = 0;
                double averagev = 0;
                long Tbias = 0;
                for (i = 0; i < seqsiz; ++i)
                    {
                    Tbias += Ts[i];
                    }
                Tbias /= seqsiz;
                for (i = index - seqsiz; i < index; ++i)
                    {
                    int id = i % seqsiz;
                    averageh += hs[id];
                    averagev += vs[id];
                    ts[id] = (double)(Ts[id] - Tbias); // This makes ts[id] a small number around 0
                    averageTime += ts[id]; // later used to slightly adjust the ts[id] so that their sum becomes zero.
                    }
                averageh = averageh / period;
                averagev = averagev / period;
                averageTime = averageTime / period;
                for (i = index - seqsiz; i < index; ++i)
                    {
                    int id = i % seqsiz;
                    double t = ts[id] - averageTime;
                    double t2 = t * t;
                    double t3 = t2 * t;
                    double h = hs[id] - averageh;
                    double v = vs[id] - averagev;
                    St += t;
                    St2 += t2;
                    St3 += t2 * t;
                    St4 += t2 * t2;
                    St5 += t2 * t3;
                    St6 += t3 * t3;
                    Sh += h;
                    Sth += t * h;
                    St2h += t2 * h;
                    St3h += t3 * h;
                    Sv += v;
                    Stv += t * v;
                    St2v += t2 * v;
                    St3v += t3 * v;
                    }
                switch (what)
                    {
                    case evelo:
                        assignToPers->vx = velocity(St, St2, Sh, Sth, period);
                        assignToPers->vy = velocity(St, St2, Sv, Stv, period);
                        break;
                    case eacce:
                        assignToPers->ax = acceleration(St, St2, St3, St4, Sh, Sth, St2h, period);
                        assignToPers->ay = acceleration(St, St2, St3, St4, Sv, Stv, St2v, period);
                        break;
                    case ejerk:
                        assignToPers->jx = jerk(St2, St2h, St3, St3h, St4, St5, St6, Sth, period);
                        assignToPers->jy = jerk(St2, St2v, St3, St3v, St4, St5, St6, Stv, period);
                        break;
                    }
                }
            else
                {
                //Not enough data collected to fill a period.
                }
            }
        }
    }

line* readInput(char* name)
    {
    FILE* fp;
    printf("name %s\n", name);
    fp = fopen(name, "r");
    if (fp)
        {
        int lines = 1;
        int kar;
        int columns = 0;
        int mincols = 99999999;
        int maxcols = 0;
        int cols = 0;
        int maxlen = 0;
        int len = 0;
        while ((kar = fgetc(fp)) != EOF)
            {
            ++len;
            if (kar == '\n')
                {
                ++lines;
                if (cols > maxcols)
                    maxcols = cols;
                if (cols < mincols)
                    mincols = cols;
                if (len > maxlen)
                    maxlen = len;
                len = 0;
                cols = 0;
                }
            else if (kar == '\t')
                {
                ++cols;
                }
            }
        if (mincols == maxcols)
            {
            columns = maxcols;
            if (columns != 3 * persons)
                {
                printf("Persons %d, columns %d\n", persons, columns);
                exit(2);
                }
            }
        else
            {
            printf("Different number of columns %d and %d\n", mincols, maxcols);
            exit(1);
            }
        rewind(fp);
        // 140	94.7639	90.7244	0.357086	590.125	110.799	0.275333
        line* Lines = new line[lines];
        char* buffer = new char[maxlen + 2];
        int L = 0;
        for (L = 0; L < lines - 1; ++L)
            {
            fgets(buffer, maxlen + 1, fp);
            int len = strlen(buffer);
            if (len > 0)
                buffer[len - 1] = '\t'; // overwrite \n
            char* tab = strchr(buffer, '\t');
            *tab++ = '\0';
            Lines[L].frame = strtol(buffer, 0L, 10);
            for (int i = 0; i < persons; ++i)
                {
                char* nexttab = strchr(tab, '\t');
                *nexttab = '\0';
                Lines[L].P[i].x = strtod(tab, 0L);
                tab = nexttab + 1;

                nexttab = strchr(tab, '\t');
                *nexttab = '\0';
                Lines[L].P[i].y = strtod(tab, 0L);
                tab = nexttab + 1;

                nexttab = strchr(tab, '\t');
                if (nexttab)
                    *nexttab = '\0';
                Lines[L].P[i].p = strtod(tab, 0L);
                if (nexttab)
                    tab = nexttab + 1;
                }
            }
        Lines[L].frame = 0;
        for (int i = 0; i < persons; ++i)
            {
            Lines[L].P[i].x = 0.0;
            Lines[L].P[i].y = 0.0;
            Lines[L].P[i].p = 0.0;
            }
        fclose(fp);
        return Lines;
        }
    else
        return 0;
    }

int clock(double sagittaH, double sagittaV)
    {
    double theta = atan2((double)-sagittaV, (double)sagittaH) * 180.0 / M_PI;
    theta = -theta;
    if (theta < 0.0)
        theta += 360;
    theta += 105; // half past eleven
    int direction = (int)theta;
    direction %= 360;
    direction /= 30; // 0 .. 11
    if (direction == 0)
        direction = 12;
    return direction;
    }

double rho(double sagittaH, double sagittaV)
    {
    double rho = sqrt((double)(sagittaV * sagittaV) + (double)(sagittaH * sagittaH));
    return rho;
    }

void printPerson(FILE* fp, struct person* Person, char* name)
    {
    fprintf(fp,
        "%lf\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf",
        Person->x, Person->y, Person->p,
        rho(Person->vx, Person->vy), clock(Person->vx, Person->vy), Person->vx, Person->vy,
        rho(Person->ax, Person->ay), clock(Person->ax, Person->ay), Person->ax, Person->ay,
        rho(Person->jx, Person->jy), clock(Person->jx, Person->jy), Person->jx, Person->jy
    );
    // F2:x-pos	F2:y-pos	F2:weight	F2:velocity-r	F2:velocity-clock	F2:velocity-x	F2:velocity-y	F2:acceleration-r	F2:acceleration-clock	F2:acceleration-x	F2:acceleration-y	F2:jerk-r	F2:jerk-clock	F2:jerk-x	F2:jerk-y
    }

void printPersonHead(FILE* fp, char* name)
    {
    fprintf(fp
        , "%s:x-pos\t%s:y-pos\t%s:weight\t%s:velocity-r\t%s:velocity-clock\t%s:velocity-x\t%s:velocity-y\t%s:acceleration-r\t%s:acceleration-clock\t%s:acceleration-x\t%s:acceleration-y\t%s:jerk-r\t%s:jerk-clock\t%s:jerk-x\t%s:jerk-y"
        , name, name, name, name, name, name, name, name, name, name, name, name, name, name, name
    );
    }

//static char names[2][3];
static char** names;

void print(struct line* Lines, const char* out)
    {
    FILE* fp = fopen(out, "wb");
    if (fp)
        {
        fprintf(fp, "frame\t");
        for (int p = 0; p < persons; ++p)
            {
            printPersonHead(fp, names[p]);
            if (p + 1 < persons)
                fputc('\t', fp);
            else
                fputc('\n', fp);
            }
        struct line* current = Lines;
        for (current = Lines; current->frame; ++current)
            {
            fprintf(fp, "%d\t", (current->frame) << 2);
            for (int p = 0; p < persons; ++p)
                {
                printPerson(fp, current->P + p, names[p]);
                if (p + 1 < persons)
                    fputc('\t', fp);
                else
                    fputc('\n', fp);
                }
            }
        fclose(fp);
        }
    }

int main(int argc, char** argv)
    {
    if (argc < 8)
        {
        //       0      1       2          3                    4                  5              6          7       8     9..
        printf("jerk <input> <output> <velocity period> <acceleration period> <jerk period> <margin type> <point> [<point> ...] \n");
        return -1;
        }
    /*
    names[0][0] = argv[1][0];
    names[0][1] = argv[1][1];
    names[0][2] = 0;
    names[1][0] = argv[1][3];
    names[1][1] = argv[1][4];
    names[1][2] = 0;
    */
    persons = argc - 7;
    names = new char* [persons];
    for (int i = 0; i < persons; i++)
        names[i] = argv[7 + i];
    struct line* Lines = readInput(argv[1]);
    if (Lines)
        {
        unsigned long velow = strtoul(argv[3], NULL, 10);
        unsigned long accew = strtoul(argv[4], NULL, 10);
        unsigned long jerkw = strtoul(argv[5], NULL, 10);
        int pers = 0;
        for (pers = 0; pers < persons; ++pers)
            {
            doPerson(pers, evelo, argv[3], argv[6], Lines);
            doPerson(pers, eacce, argv[4], argv[6], Lines);
            doPerson(pers, ejerk, argv[5], argv[6], Lines);
            }
        print(Lines, argv[2]);
        free(Lines);
        }
    return 0;
    }
